using Plotly, SmithNormalForm, LinearAlgebra, OrderedCollections, SparseArrays, Combinatorics
using Revise
using Zipper
using DataFrames
using SparseArrays

function sitefock(site::Offset; flavorcount::Integer = 1)::FockSpace{Offset}
    basis::Offset = site|>basispoint
    offset::Offset = site - basis
    return FockSpace((Mode(:offset => offset, :pos => basis, :flavor => f) for f in 1:flavorcount), reflected=site)
end

regionfock(region::Region; flavorcount::Integer = 1)::FockSpace{Region} = (
    fockspaceunion(sitefock(r, flavorcount=flavorcount) for r in region)|>FockSpace{Region})

function circularregionmodes(correlations::FockMap, origin::Offset, physicalmodes::Subset{Mode}, radius::Number)::Subset{Mode}
    currentspace::RealSpace = correlations.outspace |> getcrystal |> getspace |> orthospace
    physicalnorm = m -> lineartransform(currentspace, m |> getpos) |> norm
    return filter(p -> physicalnorm(p - origin) < radius, physicalmodes)
end

function realphyspos(mode::Mode)
    currentspace::RealSpace = mode|> getspace |> euclidean
    return(lineartransform(currentspace, mode |> getpos))
end

function _findlocalspstates(;
    statecorrelations::FockMap, regionfock::FockSpace,
    symmetry::AffineTransform = identitytransform(statecorrelations |> getcrystal |> dimension),
    spectrumextractpredicate::Function = v -> v < 1e-2,
    linearindependencethreshold::Real = 5e-2,
    degeneracythreshold::Real = 1e-7)

    function lineardependencefilter(spstate::FockMap)::Bool
        pseudoidentity::FockMap = (spstate' * spstate)
        mineigenvalue = minimum(v for (_, v) in pseudoidentity |> eigvalsh)
        return mineigenvalue > linearindependencethreshold
    end

    localcorrelations::FockMap = regioncorrelations(statecorrelations, regionfock)
    localspectrum::EigenSpectrum = eigspech(localcorrelations, groupingthreshold=degeneracythreshold)
    groupeigenvalues::Base.Generator = (
        subset => (localspectrum |> geteigenvalues)[subset |> first]
        for subset in localspectrum |> geteigenvectors |> getinspace |> sparsegrouping(:eigenindex) |> rep)
    selectedgroups = Iterators.filter(p -> p.second |> spectrumextractpredicate, groupeigenvalues)

    selectedisometries = ((localspectrum |> geteigenvectors)[:, group.first |> FockSpace] for group in selectedgroups)
    orthogonalspstates = Iterators.filter(lineardependencefilter, selectedisometries)
    symmetricspstates = (state * *(state, symmetry) for state in orthogonalspstates)
    spstates = (state * spatialmap(state)' for state in symmetricspstates)

    return (state |> getinspace |> dimension => state for state in spstates)
end

# Need further improvement
function groupmodesbydist(;
    region::Subset{Point{RealSpace}},
    regionfock::FockSpace,
    center::Point,
    samedistancethreshold::Int = 4)

    visualspace = region |> getspace |> euclidean
    physicalnorm = m -> lineartransform(visualspace, m |> getpos) |> norm
    distancewithmode = sort([(physicalnorm(mode-center),mode) for mode in regionfock], by = first, rev = true)
    df = DataFrame()
    df.distance = [round(dist; digits=samedistancethreshold) for (dist,_) in distancewithmode]
    df.mode = [mode for (_,mode) in distancewithmode]
    grouped_df = groupby(df, :distance)
    store = Dict()
    for (ind,group) in enumerate(grouped_df)
        store[ind] = []
        for (distance,mode) in zip(group.distance,group.mode)
            push!(store[ind],(distance,mode))
        end
    end
    return store
end

function _spatialmap(fockmap::FockMap)::FockMap
    function _spatialinmode(colmap::FockMap, ind::Integer)
        inmode::Mode = colmap |> getinspace |> first
        absmap::FockMap = colmap |> abs 
        modecenter::Offset = sort(absmap |> Zipper.columnspec, by=p->p.second |> real) |> last |> first |> getpos
        basis::Offset = modecenter |> basispoint
        offset::Offset = modecenter - basis
        return inmode |> setattr(:offset => offset) |> setattr(:pos => basis) 
    end

    spatialinspace::FockSpace{Region} = FockSpace{Region}( _spatialinmode(fockmap[:, m],i) for (i,m) in fockmap |> getinspace |> enumerate)
    return idmap(spatialinspace, fockmap |> getinspace)
end

function _localwannierization(localbasis::FockMap, localseeds::FockMap, svdorthothreshold::Number = 1e-1)::FockMap
    U, Σ, Vt = svd(localbasis'*localseeds)
    minsvdvalue::Number = minimum(v for (_, v) in Σ)
    println("min svdvalue", minsvdvalue)
    precarioussvdvalues::Vector = []
    if minsvdvalue < svdorthothreshold
        push!(precarioussvdvalues, minsvdvalue)
    end
    unitary::FockMap = U * Vt
    wannierizedbasis = localbasis*unitary
    wannierizedbasis = wannierizedbasis*_spatialmap(wannierizedbasis)'
    return wannierizedbasis
end

function _crystalisometries(; localisometry::FockMap, crystalfock::FockSpace{Crystal},
    addinspacemomentuminfo::Bool = false)::Dict{Momentum, FockMap}

    crystal::Crystal = getcrystal(crystalfock)
    # fouriermap::FockMap = fourier(crystalfock, localisometry.outspace) / (crystal |> vol |> sqrt) seems do not need the denominator
    fouriermap::FockMap = fourier(crystalfock, localisometry.outspace) 
    momentumfouriers::Base.Generator = rowsubmaps(fouriermap)
    bz::Subset{Momentum} = brillouinzone(crystal)


    function preprocesslocalisometry(k::Momentum)::FockMap
        if !addinspacemomentuminfo
            return localisometry
        end
        inspace::FockSpace = localisometry.inspace |> orderedmodes |> setattr(:offset => k) |> FockSpace
        return FockMap(localisometry, inspace=inspace, performpermute=false)
    end

    return Dict(k => kfourier * preprocesslocalisometry(k) for (k, kfourier) in zip(bz, momentumfouriers))
end

function constructfilledandemptyseed(modewithdistlist,symmetry::AffineTransform)
    firstpos = modewithdistlist[1][2]|>getpos
    emptyseedspos = [firstpos,symmetry^2*firstpos,symmetry^4*firstpos]
    filledseedspos = [symmetry*firstpos,symmetry^3*firstpos,symmetry^5*firstpos]
    return emptyseedspos, filledseedspos
end

triangular = RealSpace([sqrt(3)/2 -1/2; 0. 1.]')
kspace = convert(MomentumSpace, triangular)

pa = [1/3, 2/3] ∈ triangular
pb = [2/3, 1/3] ∈ triangular
pc = (pa + pb) / 2
spatialsnappingcalibration((pa, pb, pc))

c6 = pointgrouptransform([cos(π/3) -sin(π/3); sin(π/3) cos(π/3)])

unitcell = Subset(pa, pb)
crystal = Crystal(unitcell, [24, 24])

reciprocalhashcalibration(crystal.sizes)

modes::Subset{Mode} = quantize(:pos, unitcell, 1)
m0, m1 = members(modes)

tₙ = -1 + 0im
tₕ = 0.1im

nearestneighbor = [
    (m0, m1) => tₙ,
    (m0, m1 |> setattr(:offset => [-1, 0] ∈ triangular)) => tₙ,
    (m0, m1 |> setattr(:offset => [0, 1] ∈ triangular)) => tₙ]

haldane = [
    (m0, setattr(m0, :offset => Point([1, 1], triangular))) => tₕ,
    (m0, setattr(m0, :offset => Point([-1, 0], triangular))) => tₕ,
    (m0, setattr(m0, :offset => Point([0, -1], triangular))) => tₕ,
    (m1, setattr(m1, :offset => Point([1, 1], triangular))) => -tₕ,
    (m1, setattr(m1, :offset => Point([-1, 0], triangular))) => -tₕ,
    (m1, setattr(m1, :offset => Point([0, -1], triangular))) => -tₕ]

bonds::FockMap = bondmap([nearestneighbor..., haldane...])

energyspectrum = computeenergyspectrum(bonds, crystal=crystal)
energyspectrum |> visualize

groundstates::CrystalSpectrum = groundstatespectrum(energyspectrum, perunitcellfillings=1)
groundstates |> visualize

groundstateprojector = groundstates |> crystalprojector
C = idmap(groundstateprojector.outspace) - groundstateprojector
correlations = C
CC = CrystalFockMap(C)
espec = correlations|>eigspec

crystalfock = correlations.outspace

scale = Scale([2 0; 0 2], crystalfock |> getcrystal |> getspace)
@time blockresult = blocking(:action => scale, :correlations => correlations, :crystal => (crystalfock |> getcrystal))

blockedcrystal::Crystal = blockresult[:crystal]
blockedcorrelations::FockMap = blockresult[:correlations]
blocker = blockresult[:transformer]

crystalpoints::Subset{Offset} = latticepoints(blockedcrystal)
samplepoints::Subset{Offset} = crystalpoints + c6^2 * crystalpoints + c6^4 * crystalpoints
blockedmodes::Subset{Mode} = quantize(:pos, blockedcrystal.unitcell, 1)
physicalmodes::Subset{Mode} = spanoffset(blockedmodes, samplepoints)

localmodes = circularregionmodes(blockedcorrelations, [0, 0] ∈ (blockedmodes |> getspace) , physicalmodes, 1)
localfock = FockSpace{Region}(localmodes)
visualize(getregion(localfock),visualspace=euclidean(RealSpace, 2))
visualize(regioncorrelations(blockedcorrelations, localfock) |> eigspec)

blockedcrystal = blockedcorrelations.outspace |> getcrystal
scaleforgmera = Scale([3 0; 0 3], blockedcrystal |> getspace)
@time scaleresult = blocking(:action => scaleforgmera , :correlations => blockedcorrelations, :crystal => blockedcrystal)

function onestepgmera(correlations::FockMap,localcenter1input,localcenter2input,localcenter3input)
    crystal::Crystal = getcrystal(correlations.inspace)

    crystalpoints::Subset{Offset} = latticepoints(crystal)
    samplepoints::Subset{Offset} = crystalpoints + c6^2 * crystalpoints + c6^4 * crystalpoints
    modes::Subset{Mode} = quantize(:pos, crystal.unitcell, 1)
    physicalmodes::Subset{Mode} = spanoffset(modes, samplepoints)

    localcenter1::Offset = localcenter1input ∈ (modes |> getspace)
    localcenter2::Offset = localcenter2input ∈ (modes |> getspace)
    localcenter3::Offset = localcenter3input ∈ (modes |> getspace)


    function fulllocalwannierization(correlations::FockMap,localcenter::Offset,radius::Real=1/3,selectionstragedy::Function=modeselectionbycount(3))
        localmodes = circularregionmodes(correlations, localcenter , physicalmodes, radius)
        localfock::FockSpace = FockSpace{Region}(localmodes)
        localregion::Subset{Offset} = Subset(m |> getpos for m in localmodes)
        modebydist = groupmodesbydist(region = localregion,regionfock = localfock,center = localcenter) 
        # group the eigenvectors to isometries with filled, empty ,courier
        localiso = localisometries(correlations, localfock, selectionstrategy=selectionstragedy)
        # finding seeds for local frozen (6 modes at the center)
        frozenseedingfocks = (FockSpace{Region}(circularregionmodes(correlations, modewithdist[2]|> getpos, physicalmodes, 0.01)) for modewithdist in modebydist[3])
        frozenseeds = sum([reduce(+ ,FockMap(v,inspace = v |> getinspace |> orderedmodes |> FockSpace,performpermute = false) for (n, (_ , v)) in enumerate(_findlocalspstates(statecorrelations = correlations, regionfock = frozenseedingfock, 
        spectrumextractpredicate = v -> true, symmetry = identitytransform(2)))) for frozenseedingfock in frozenseedingfocks])
        extendedfrozenseeds = idmap(localfock, localfock)[:,frozenseeds |> getoutspace] * frozenseeds
        wannierizedfrozens = _localwannierization(localiso[:filled]+localiso[:empty], extendedfrozenseeds)

        emptyseedspos, filledseedspos = constructfilledandemptyseed(modebydist[3],c6|> recenter(localcenter))

        emptyseedingfocks = (FockSpace{Region}(circularregionmodes(correlations, pos, physicalmodes, 0.01)) for pos in emptyseedspos)
        emptyseeds = sum([reduce(+ ,FockMap(v,inspace = v |> getinspace |> orderedmodes |> FockSpace,performpermute = false) for (n, (_ , v)) in enumerate(_findlocalspstates(statecorrelations = correlations, regionfock = emptyseedingfock, 
        spectrumextractpredicate = v -> true, symmetry = identitytransform(2)))) for emptyseedingfock in emptyseedingfocks])
        extendedemptyseeds = idmap(localfock, localfock)[:,emptyseeds |> getoutspace] * emptyseeds
        wannierizedemptys = _localwannierization(localiso[:empty], extendedemptyseeds)

        filledseedingfocks = (FockSpace{Region}(circularregionmodes(correlations, pos, physicalmodes, 0.01)) for pos in filledseedspos)
        filledseeds = sum([reduce(+ ,FockMap(v,inspace = v |> getinspace |> orderedmodes |> FockSpace,performpermute = false) for (n, (_ , v)) in enumerate(_findlocalspstates(statecorrelations = correlations, regionfock = filledseedingfock, 
        spectrumextractpredicate = v -> true, symmetry = identitytransform(2)))) for filledseedingfock in filledseedingfocks])
        extendedfilledseeds = idmap(localfock, localfock)[:,filledseeds |> getoutspace] * filledseeds
        wannierizedfilleds = _localwannierization(localiso[:filled], extendedfilledseeds)

        
        return wannierizedemptys, wannierizedfilleds, wannierizedfrozens, frozenseeds, filledseeds, emptyseeds
    end

    wannierizedemptys1, wannierizedfilleds1, wannierizedfrozens1, frozenseeds1, filledseeds1, emptyseeds1 = fulllocalwannierization(correlations, localcenter1)
    wannierizedemptys2, wannierizedfilleds2, wannierizedfrozens2, frozenseeds2, filledseeds2, emptyseeds2 = fulllocalwannierization(correlations, localcenter2)
    wannierizedemptys3, wannierizedfilleds3, wannierizedfrozens3, frozenseeds3, filledseeds3, emptyseeds3 = fulllocalwannierization(correlations, localcenter3)

    localwannierizedfrozens = wannierizedfrozens1 + wannierizedfrozens2 + wannierizedfrozens3
    localwannierizedemptys = wannierizedemptys1 + wannierizedemptys2 + wannierizedemptys3 
    localwannierizedfilleds = wannierizedfilleds1 + wannierizedfilleds2 + wannierizedfilleds3

    frozenseeds = frozenseeds1 + frozenseeds2 + frozenseeds3
    emptyseeds = emptyseeds1 + emptyseeds2 + emptyseeds3
    filledseeds = filledseeds1 + filledseeds2 + filledseeds3

    function globalwannierfunction(seeds, localisometry::FockMap)
        wanniercrystalisos = _crystalisometries(localisometry=localisometry, crystalfock=correlations.outspace, addinspacemomentuminfo = true)
        wannierunitcell::Subset{Offset} = Subset((mode |> getattr(:pos))  for mode in seeds.inspace |> orderedmodes)
        
        wanniercrystall::Crystal = Crystal(wannierunitcell, (correlations.outspace |> getcrystal).sizes)
        # visualize(latticepoints(wanniercrystall))
        
        wannierisometry::FockMap = directsum(wanniercrystaliso for (k,wanniercrystaliso) in wanniercrystalisos)
        inspace::FockSpace = FockSpace(wannierisometry.inspace, reflected=wanniercrystall)
        outspace::FockSpace = FockSpace(wannierisometry.outspace, reflected=scaledcorrelations.outspace)
        globalwannierizedfunction = FockMap(wannierisometry, inspace=inspace, outspace=outspace, performpermute=false)

        return globalwannierizedfunction
    end

    wannierfrozenisometry = globalwannierfunction(frozenseeds,localwannierizedfrozens)
    wannieremptyisometry = globalwannierfunction(emptyseeds,localwannierizedemptys)
    wannierfilledisometry = globalwannierfunction(filledseeds,localwannierizedfilleds)

    frozencorrelations = wannierfrozenisometry' * correlations * wannierfrozenisometry
    emptycorrelations = wannieremptyisometry' * correlations * wannieremptyisometry
    filledcorrelations = wannierfilledisometry' * correlations * wannierfilledisometry

    return Dict(
        :emptyisometry => wannieremptyisometry,
        :filledisometry => wannierfilledisometry,
        :frozenisometry => wannierfrozenisometry,
        :filledcorrelations => filledcorrelations, 
        :emptycorrelations => emptycorrelations,
        :frozencorrelations => frozencorrelations)

end

localcenter1input = [1/3, 0] 
localcenter2input = [0, 1/3] 
localcenter3input = [-1/3, -1/3]

onestepgmera(scaleresult[:correlations],localcenter1input,localcenter2input,localcenter3input)



scaledcrystal = scaleresult[:crystal]
scaledcorrelations = scaleresult[:correlations]
scaler = scaleresult[:transformer]

scaledcrystal = scaledcorrelations.inspace |> getcrystal

refcrystalpoints::Subset{Offset} = latticepoints(scaledcrystal)
refsamplepoints::Subset{Offset} = refcrystalpoints + c6^2 * refcrystalpoints + c6^4 * refcrystalpoints
refblockedmodes::Subset{Mode} = quantize(:pos, scaledcrystal.unitcell, 1)
refphysicalmodes::Subset{Mode} = spanoffset(refblockedmodes, refsamplepoints)

localcenter1::Offset = [0, 0] ∈ (refblockedmodes |> getspace)
localmodes1 = circularregionmodes(scaledcorrelations, localcenter1 , refphysicalmodes, 1/3)
localfock1::FockSpace = FockSpace{Region}(localmodes1)
localregion1::Subset{Offset} = Subset(m |> getpos for m in localmodes1)
modebydist1 = groupmodesbydist(region = localregion1,regionfock = localfock1,center = localcenter1) 
# group the eigenvectors to isometries with filled, empty ,courier
localiso1 = localisometries(scaledcorrelations, localfock1, selectionstrategy=modeselectionbycount(3))
# finding seeds for local frozen (6 modes at the center)
frozenseedingfocks1 = (FockSpace{Region}(circularregionmodes(scaledcorrelations, modewithdist[2]|> getpos, refphysicalmodes, 0.01)) for modewithdist in modebydist1[3])
frozenseeds1 = sum([reduce(+ ,FockMap(v,inspace = v |> getinspace |> orderedmodes |> FockSpace,performpermute = false) for (n, (_ , v)) in enumerate(_findlocalspstates(statecorrelations = scaledcorrelations, regionfock = frozenseedingfock, 
spectrumextractpredicate = v -> true, symmetry = identitytransform(2)))) for frozenseedingfock in frozenseedingfocks1])
extendedfrozenseeds1 = idmap(localfock1, localfock1)[:,frozenseeds1 |> getoutspace] * frozenseeds1
wannierizedfrozens1 = _localwannierization(localiso1[:filled]+localiso1[:empty], extendedfrozenseeds1)

emptyseedspos1, filledseedspos1 = constructfilledandemptyseed(modebydist1[3],c6|> recenter(localcenter1))

emptyseedingfocks1 = (FockSpace{Region}(circularregionmodes(scaledcorrelations, pos, refphysicalmodes, 0.01)) for pos in emptyseedspos1)
emptyseeds1 = sum([reduce(+ ,FockMap(v,inspace = v |> getinspace |> orderedmodes |> FockSpace,performpermute = false) for (n, (_ , v)) in enumerate(_findlocalspstates(statecorrelations = scaledcorrelations, regionfock = emptyseedingfock, 
spectrumextractpredicate = v -> true, symmetry = identitytransform(2)))) for emptyseedingfock in emptyseedingfocks1])
extendedemptyseeds1 = idmap(localfock1, localfock1)[:,emptyseeds1 |> getoutspace] * emptyseeds1
wannierizedempty1 = _localwannierization(localiso1[:empty], extendedemptyseeds1)

filledseedingfocks1 = (FockSpace{Region}(circularregionmodes(scaledcorrelations, pos, refphysicalmodes, 0.01)) for pos in filledseedspos1)
filledseeds1 = sum([reduce(+ ,FockMap(v,inspace = v |> getinspace |> orderedmodes |> FockSpace,performpermute = false) for (n, (_ , v)) in enumerate(_findlocalspstates(statecorrelations = scaledcorrelations, regionfock = filledseedingfock, 
spectrumextractpredicate = v -> true, symmetry = identitytransform(2)))) for filledseedingfock in filledseedingfocks1])
extendedfilledseeds1 = idmap(localfock1, localfock1)[:,filledseeds1 |> getoutspace] * filledseeds1
wannierizedfilled1 = _localwannierization(localiso1[:filled], extendedfilledseeds1)

visualize(wannierizedfrozens1'*regioncorrelations(scaledcorrelations, localfock1)*wannierizedfrozens1|>eigspec)


localcenter2::Offset = [1/3, 2/3] ∈ (refblockedmodes |> getspace)
localmodes2 = circularregionmodes(scaledcorrelations, localcenter2 , refphysicalmodes, 1/3)
localfock2::FockSpace = FockSpace{Region}(localmodes2)
localregion2::Subset{Offset} = Subset(m |> getpos for m in localmodes2)
modebydist2 = groupmodesbydist(region = localregion2,regionfock = localfock2,center = localcenter2)  
# group the eigenvectors to isometries with filled, empty ,courier
localiso2 = localisometries(scaledcorrelations, localfock2, selectionstrategy=modeselectionbycount(3))
# finding seeds for local frozen (6 modes at the center)
frozenseedingfocks2 = (FockSpace{Region}(circularregionmodes(scaledcorrelations, modewithdist[2]|> getpos, refphysicalmodes, 0.01)) for modewithdist in modebydist2[3])
frozenseeds2 = sum([reduce(+ ,FockMap(v,inspace = v |> getinspace |> orderedmodes |> FockSpace,performpermute = false) for (n, (_ , v)) in enumerate(_findlocalspstates(statecorrelations = scaledcorrelations, regionfock = frozenseedingfock, 
spectrumextractpredicate = v -> true, symmetry = identitytransform(2)))) for frozenseedingfock in frozenseedingfocks2])
extendedfrozenseeds2 = idmap(localfock2, localfock2)[:,frozenseeds2 |> getoutspace] * frozenseeds2
wannierizedfrozens2 = _localwannierization(localiso2[:filled]+localiso2[:empty], extendedfrozenseeds2)

emptyseedspos2, filledseedspos2 = constructfilledandemptyseed(modebydist2[3],c6|> recenter(localcenter2))

emptyseedingfocks2 = (FockSpace{Region}(circularregionmodes(scaledcorrelations, pos, refphysicalmodes, 0.01)) for pos in emptyseedspos2)
emptyseeds2 = sum([reduce(+ ,FockMap(v,inspace = v |> getinspace |> orderedmodes |> FockSpace,performpermute = false) for (n, (_ , v)) in enumerate(_findlocalspstates(statecorrelations = scaledcorrelations, regionfock = emptyseedingfock, 
spectrumextractpredicate = v -> true, symmetry = identitytransform(2)))) for emptyseedingfock in emptyseedingfocks2])
extendedemptyseeds2 = idmap(localfock2, localfock2)[:,emptyseeds2 |> getoutspace] * emptyseeds2
wannierizedempty2 = _localwannierization(localiso2[:empty], extendedemptyseeds2)

filledseedingfocks2 = (FockSpace{Region}(circularregionmodes(scaledcorrelations, pos, refphysicalmodes, 0.01)) for pos in filledseedspos2)
filledseeds2 = sum([reduce(+ ,FockMap(v,inspace = v |> getinspace |> orderedmodes |> FockSpace,performpermute = false) for (n, (_ , v)) in enumerate(_findlocalspstates(statecorrelations = scaledcorrelations, regionfock = filledseedingfock, 
spectrumextractpredicate = v -> true, symmetry = identitytransform(2)))) for filledseedingfock in filledseedingfocks2])
extendedfilledseeds2 = idmap(localfock2, localfock2)[:,filledseeds2 |> getoutspace] * filledseeds2
wannierizedfilled2 = _localwannierization(localiso2[:filled], extendedfilledseeds2)

visualize(wannierizedempty2'*regioncorrelations(scaledcorrelations, localfock2)*wannierizedempty2|>eigspec)

localcenter3::Offset = [-1/3, 1/3] ∈ (refblockedmodes |> getspace)
localmodes3 = circularregionmodes(scaledcorrelations, localcenter3 , refphysicalmodes, 1/3)
localfock3::FockSpace = FockSpace{Region}(localmodes3)
localregion3::Subset{Offset} = Subset(m |> getpos for m in localmodes3)
modebydist3 = groupmodesbydist(region = localregion3,regionfock = localfock3,center = localcenter3)
# group the eigenvectors to isometries with filled, empty ,courier
localiso3 = localisometries(scaledcorrelations, localfock3, selectionstrategy=modeselectionbycount(3))
# finding seeds for local frozen (6 modes at the center)
frozenseedingfocks3 = (FockSpace{Region}(circularregionmodes(scaledcorrelations, modewithdist[2]|> getpos, refphysicalmodes, 0.01)) for modewithdist in modebydist3[3])
frozenseeds3 = sum([reduce(+ ,FockMap(v,inspace = v |> getinspace |> orderedmodes|> FockSpace,performpermute = false) for (n, (_ , v)) in enumerate(_findlocalspstates(statecorrelations = scaledcorrelations, regionfock = frozenseedingfock, 
spectrumextractpredicate = v -> true, symmetry = identitytransform(2)))) for frozenseedingfock in frozenseedingfocks3])
extendedfrozenseeds3 = idmap(localfock3, localfock3)[:,frozenseeds3 |> getoutspace] * frozenseeds3
wannierizedfrozens3 = _localwannierization(localiso3[:filled]+localiso3[:empty], extendedfrozenseeds3)

emptyseedspos3, filledseedspos3 = constructfilledandemptyseed(modebydist3[3],c6|> recenter(localcenter3))

emptyseedingfocks3 = (FockSpace{Region}(circularregionmodes(scaledcorrelations, pos, refphysicalmodes, 0.01)) for pos in emptyseedspos3)
emptyseeds3 = sum([reduce(+ ,FockMap(v,inspace = v |> getinspace |> orderedmodes |> FockSpace,performpermute = false) for (n, (_ , v)) in enumerate(_findlocalspstates(statecorrelations = scaledcorrelations, regionfock = emptyseedingfock, 
spectrumextractpredicate = v -> true, symmetry = identitytransform(2)))) for emptyseedingfock in emptyseedingfocks3])
extendedemptyseeds3 = idmap(localfock3, localfock3)[:,emptyseeds3 |> getoutspace] * emptyseeds3
wannierizedempty3 = _localwannierization(localiso3[:empty], extendedemptyseeds3)

filledseedingfocks3 = (FockSpace{Region}(circularregionmodes(scaledcorrelations, pos, refphysicalmodes, 0.01)) for pos in filledseedspos3)
filledseeds3 = sum([reduce(+ ,FockMap(v,inspace = v |> getinspace |> orderedmodes |> FockSpace,performpermute = false) for (n, (_ , v)) in enumerate(_findlocalspstates(statecorrelations = scaledcorrelations, regionfock = filledseedingfock, 
spectrumextractpredicate = v -> true, symmetry = identitytransform(2)))) for filledseedingfock in filledseedingfocks3])
extendedfilledseeds3 = idmap(localfock3, localfock3)[:,filledseeds3 |> getoutspace] * filledseeds3
wannierizedfilled3 = _localwannierization(localiso3[:filled], extendedfilledseeds3)

visualize(wannierizedempty3'*regioncorrelations(scaledcorrelations, localfock3)*wannierizedempty3|>eigspec)

visualize(wannierizedfrozens3'*regioncorrelations(scaledcorrelations, localfock3)*wannierizedfrozens3|>eigspec)

# visualize((wannierizedfrozens3'*regioncorrelations(scaledcorrelations, localfock3)*wannierizedfrozens3) |> eigspec)

frozenseeds = frozenseeds1 + frozenseeds2 + frozenseeds3
emptyseeds = emptyseeds1 + emptyseeds2 + emptyseeds3
filledseeds = filledseeds1 + filledseeds2 + filledseeds3

wannierizedfrozens = wannierizedfrozens1 + wannierizedfrozens2 + wannierizedfrozens3
wannierizedemptys = wannierizedempty1 + wannierizedempty2 + wannierizedempty3
wannierizedfilleds = wannierizedfilled1 + wannierizedfilled2 + wannierizedfilled3


wanniercrystalisos = _crystalisometries(localisometry=wannierizedfrozens, crystalfock=scaledcorrelations.outspace, addinspacemomentuminfo = true)
wannierunitcell::Subset{Offset} = Subset((mode |> getattr(:pos))  for mode in frozenseeds.inspace |> orderedmodes)
visualize(wannierunitcell)
visualize(latticepoints(scaledcorrelations.outspace |> getcrystal))

wanniercrystall::Crystal = Crystal(wannierunitcell, (scaledcorrelations.outspace |> getcrystal).sizes)
# visualize(latticepoints(wanniercrystall))

wannierisometry::FockMap = directsum(wanniercrystaliso for (k,wanniercrystaliso) in wanniercrystalisos)
inspace::FockSpace = FockSpace(wannierisometry.inspace, reflected=wanniercrystall)
outspace::FockSpace = FockSpace(wannierisometry.outspace, reflected=scaledcorrelations.outspace)
wannierizedfrozen = FockMap(wannierisometry, inspace=inspace, outspace=outspace, performpermute=false)


frozencorrelation = (wannierizedfrozen'*scaledcorrelations*wannierizedfrozen)
visualize(frozencorrelation|> eigspec)
# visualize(frozencorrelation)
Subset(getmodes(frozencorrelation.inspace))
trialphysicalmodes
visualize(regioncorrelations(frozencorrelation,FockSpace{Region}(trialphysicalmodes))|>eigspec)

(frozencorrelation.inspace |> getcrystal).unitcell
visualize((frozencorrelation.inspace |> getcrystal).unitcell)
(frozencorrelation.inspace |> getcrystal).unitcell

visualize(regioncorrelations(frozencorrelation, FockSpace{Region}(trialblockedmodes)) |> eigspec)

trialcrystalpoints::Subset{Offset} = latticepoints(frozencorrelation.inspace |> getcrystal)
trialsamplepoints::Subset{Offset} = trialcrystalpoints + c6^3 * trialcrystalpoints 
trialblockedmodes::Subset{Mode} = quantize(:pos, (frozencorrelation.inspace |> getcrystal).unitcell, 1)
trialphysicalmodes::Subset{Mode} = spanoffset(trialblockedmodes, trialsamplepoints)
trialmodes = circularregionmodes(frozencorrelation, [2/3,1/3] ∈ (trialblockedmodes |> getspace), trialphysicalmodes, 1/3)
trialfock = FockSpace{Region}(trialmodes)
visualize(getregion(trialfock),visualspace=euclidean(RealSpace, 2))
visualize(regioncorrelations(frozencorrelation, trialfock) |> eigspec)

trialmodes2 = circularregionmodes(frozencorrelation, [1/3,1/3] ∈ (trialblockedmodes |> getspace), trialphysicalmodes, 1/3)
trialfock2 = FockSpace{Region}(trialmodes2)
visualize(getregion(trialfock2),visualspace=euclidean(RealSpace, 2))
visualize(regioncorrelations(frozencorrelation, trialfock2) |> eigspec)

wanniercrystalisos = _crystalisometries(localisometry=wannierizedemptys, crystalfock=scaledcorrelations.outspace, addinspacemomentuminfo = true)
wannierunitcell::Subset{Offset} = Subset((mode |> getattr(:pos))  for mode in emptyseeds.inspace |> orderedmodes)
visualize(wannierunitcell)
visualize(latticepoints(scaledcorrelations.outspace |> getcrystal))

wanniercrystall::Crystal = Crystal(wannierunitcell, (scaledcorrelations.outspace |> getcrystal).sizes)
# visualize(latticepoints(wanniercrystall))

wannierisometry::FockMap = directsum(wanniercrystaliso for (k,wanniercrystaliso) in wanniercrystalisos)
inspace::FockSpace = FockSpace(wannierisometry.inspace, reflected=wanniercrystall)
outspace::FockSpace = FockSpace(wannierisometry.outspace, reflected=scaledcorrelations.outspace)
wannierizedempty = FockMap(wannierisometry, inspace=inspace, outspace=outspace, performpermute=false)


emptycorrelation = (wannierizedempty'*scaledcorrelations*wannierizedempty)
visualize(emptycorrelation)

(emptycorrelation.inspace |> getcrystal).unitcell
visualize((emptycorrelation.inspace |> getcrystal).unitcell)
trialcrystalpoints::Subset{Offset} = latticepoints(emptycorrelation.inspace |> getcrystal)
trialsamplepoints::Subset{Offset} = trialcrystalpoints + c6^2 * trialcrystalpoints + c6^4 * trialcrystalpoints
trialblockedmodes::Subset{Mode} = quantize(:pos, (emptycorrelation.inspace |> getcrystal).unitcell, 1)
trialphysicalmodes::Subset{Mode} = spanoffset(trialblockedmodes, trialsamplepoints)
trialmodes = circularregionmodes(emptycorrelation, [0,0] ∈ (trialblockedmodes |> getspace), trialphysicalmodes, 1)
trialfock = FockSpace{Region}(trialmodes)
visualize(getregion(trialfock),visualspace=euclidean(RealSpace, 2))
visualize(regioncorrelations(emptycorrelation, trialfock) |> eigspec)

c6*wannierizedfrozens


modeattrs(localfock)
modeattrs(blockedcorrelations.outspace)
modeattrs(trialfock)
modeattrs(frozencorrelation.outspace)
# modeattrs(newfrozencorrelation.outspace)

# new_inspace = FockSpace(mode|>setattr(:flavor=>1) |> removeattr(:ind) |> removeattr(:dumind) for mode in frozencorrelation.inspace)
# new_crystalinspace = FockSpace(new_inspace,reflected=scaledcrystal)
# newfrozencorrelation = FockMap(frozencorrelation,inspace = new_crystalinspace,outspace = new_crystalinspace,performpermute=false)


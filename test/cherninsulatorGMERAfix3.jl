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
scaler = scaleresult[:transformer]

function firststepgmera(correlations::FockMap,localcenter1input,localcenter2input,localcenter3input)
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

        # finding seeds for local courier (18 (12+6) modes at "boundary")
        courierseedingfocks1 = (FockSpace{Region}(circularregionmodes(correlations, modewithdist[2]|> getpos, physicalmodes, 0.01)) for modewithdist in modebydist[1])
        courierseedingfocks2 = (FockSpace{Region}(circularregionmodes(correlations, modewithdist[2]|> getpos, physicalmodes, 0.01)) for modewithdist in modebydist[2])
        courierseeds1 = sum([reduce(+ ,FockMap(v,inspace = v |> getinspace |> orderedmodes |> FockSpace,performpermute = false) for (n, (_ , v)) in enumerate(_findlocalspstates(statecorrelations = correlations, regionfock = courierseedingfock, 
        spectrumextractpredicate = v -> true, symmetry = identitytransform(2)))) for courierseedingfock in courierseedingfocks1])
        courierseeds2 = sum([reduce(+ ,FockMap(v,inspace = v |> getinspace |> orderedmodes |> FockSpace,performpermute = false) for (n, (_ , v)) in enumerate(_findlocalspstates(statecorrelations = correlations, regionfock = courierseedingfock, 
        spectrumextractpredicate = v -> true, symmetry = identitytransform(2)))) for courierseedingfock in courierseedingfocks2])

        combinedcourierseeds = courierseeds1 + courierseeds2
        extendedcourierseeds = idmap(localfock, localfock)[:, combinedcourierseeds |> getoutspace] * combinedcourierseeds
        wannierizedcouriers = _localwannierization(localiso[:courier], extendedcourierseeds)
        
        return wannierizedemptys, wannierizedfilleds, wannierizedfrozens, wannierizedcouriers, frozenseeds, filledseeds, emptyseeds, combinedcourierseeds
    end

    wannierizedemptys1, wannierizedfilleds1, wannierizedfrozens1, wannierizedcouriers1,frozenseeds1, filledseeds1, emptyseeds1, combinedcourierseeds1 = fulllocalwannierization(correlations, localcenter1)
    wannierizedemptys2, wannierizedfilleds2, wannierizedfrozens2, wannierizedcouriers2,frozenseeds2, filledseeds2, emptyseeds2, combinedcourierseeds2 = fulllocalwannierization(correlations, localcenter2)
    wannierizedemptys3, wannierizedfilleds3, wannierizedfrozens3, wannierizedcouriers3,frozenseeds3, filledseeds3, emptyseeds3, combinedcourierseeds3 = fulllocalwannierization(correlations, localcenter3)

    localwannierizedfrozens = wannierizedfrozens1 + wannierizedfrozens2 + wannierizedfrozens3
    localwannierizedemptys = wannierizedemptys1 + wannierizedemptys2 + wannierizedemptys3 
    localwannierizedfilleds = wannierizedfilleds1 + wannierizedfilleds2 + wannierizedfilleds3
    localwannierizedcouriers =  wannierizedcouriers1 + wannierizedcouriers2 + wannierizedcouriers3

    frozenseeds = frozenseeds1 + frozenseeds2 + frozenseeds3
    emptyseeds = emptyseeds1 + emptyseeds2 + emptyseeds3
    filledseeds = filledseeds1 + filledseeds2 + filledseeds3
    combinedcourierseeds = combinedcourierseeds1 + combinedcourierseeds2 + combinedcourierseeds3

    function globalwannierfunction(seeds, localisometry::FockMap)
        wanniercrystalisos = _crystalisometries(localisometry=localisometry, crystalfock=correlations.outspace, addinspacemomentuminfo = true)
        wannierunitcell::Subset{Offset} = Subset((mode |> getattr(:pos)) + (mode |> getattr(:offset)) for mode in seeds.inspace |> orderedmodes)
        
        wanniercrystall::Crystal = Crystal(wannierunitcell, (correlations.outspace |> getcrystal).sizes)
        # visualize(latticepoints(wanniercrystall))
        
        wannierisometry::FockMap = directsum(wanniercrystaliso for (k,wanniercrystaliso) in wanniercrystalisos)
        inspace::FockSpace = FockSpace(wannierisometry.inspace, reflected=wanniercrystall)
        outspace::FockSpace = FockSpace(wannierisometry.outspace, reflected=correlations.outspace)
        globalwannierizedfunction = FockMap(wannierisometry, inspace=inspace, outspace=outspace, performpermute=false)

        return globalwannierizedfunction
    end

    wannierfrozenisometry = globalwannierfunction(frozenseeds,localwannierizedfrozens)
    wannieremptyisometry = globalwannierfunction(emptyseeds,localwannierizedemptys)
    wannierfilledisometry = globalwannierfunction(filledseeds,localwannierizedfilleds)
    wanniercourierisometry = globalwannierfunction(combinedcourierseeds, localwannierizedcouriers)

    frozencorrelations = wannierfrozenisometry' * correlations * wannierfrozenisometry
    emptycorrelations = wannieremptyisometry' * correlations * wannieremptyisometry
    filledcorrelations = wannierfilledisometry' * correlations * wannierfilledisometry
    couriercorrelations = wanniercourierisometry' * correlations * wanniercourierisometry

    return Dict(
        :emptyisometry => wannieremptyisometry,
        :filledisometry => wannierfilledisometry,
        :frozenisometry => wannierfrozenisometry,
        :courierisometry => wanniercourierisometry,
        :filledcorrelations => filledcorrelations, 
        :emptycorrelations => emptycorrelations,
        :frozencorrelations => frozencorrelations,
        :couriercorrelations => couriercorrelations)

end


localcenter1input = [1/3, 0] 
localcenter2input = [0, 1/3] 
localcenter3input = [-1/3, -1/3]

firstgmerarea = firststepgmera(scaleresult[:correlations],localcenter1input,localcenter2input,localcenter3input)

# visualize(Subset((mode|> getattr(:offset))+(mode|> getattr(:pos))  for mode in firstgmerarea[:combinedcourierseeds].inspace |> orderedmodes))
# for mode in firstgmerarea[:combinedcourierseeds].inspace |> orderedmodes
#     println((mode|> getattr(:offset))+(mode|> getattr(:pos)))
# end

refcorrelations = firstgmerarea[:couriercorrelations]
visualize((refcorrelations.inspace |> getcrystal).unitcell)

refcrystalpoints::Subset{Offset} = latticepoints(refcorrelations.inspace |> getcrystal)
refsamplepoints::Subset{Offset} = refcrystalpoints + c6^2 * refcrystalpoints + c6^4 * refcrystalpoints
refblockedmodes::Subset{Mode} = quantize(:pos, (refcorrelations.inspace |> getcrystal).unitcell, 1)
refphysicalmodes::Subset{Mode} = spanoffset(refblockedmodes, refsamplepoints)
# visualize(getregion(FockSpace{Region}(refphysicalmodes)),visualspace=euclidean(RealSpace, 2))
trialmodes = circularregionmodes(refcorrelations, [1/3,0] ∈ (refblockedmodes |> getspace), refphysicalmodes, 1/3)
trialfock = FockSpace{Region}(trialmodes)

visualize(getregion(trialfock),visualspace=euclidean(RealSpace, 2))
visualize(regioncorrelations(refcorrelations, trialfock) |> eigspec)

# refphysicalmodes::Subset{Mode} = spanoffset(refblockedmodes, refcrystalpoints)
# visualize(getregion(FockSpace{Region}(refphysicalmodes)),visualspace=euclidean(RealSpace, 2))
trialmodes1 = circularregionmodes(refcorrelations, [1/3,0] ∈ (refblockedmodes |> getspace), refblockedmodes, 1/3)
trialfock1 = FockSpace{Region}(trialmodes1)
trialmodes2 = circularregionmodes(refcorrelations, [1/3,1] ∈ (refblockedmodes |> getspace), refblockedmodes, 1/3)
trialfock2 = FockSpace{Region}(trialmodes2)

visualize(regioncorrelations(refcorrelations, trialfock1+trialfock2) |> eigspec)
# visualize(getregion(trialfock2),visualspace=euclidean(RealSpace, 2))




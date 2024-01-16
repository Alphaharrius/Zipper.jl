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

latticevec = [sqrt(3)/2 -1/2; 0. 1.]'
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

localmodes = circularregionmodes(blockedcorrelations, [1, 0] ∈ (blockedmodes |> getspace) , physicalmodes, 1)
localfock = FockSpace{Region}(localmodes)
visualize(getregion(localfock),visualspace=euclidean(RealSpace, 2))
visualize(regioncorrelations(blockedcorrelations, localfock) |> eigspec)


blockedcrystal = blockedcorrelations.outspace |> getcrystal
scaleforgmera = Scale([3 0; 0 3], blockedcrystal |> getspace)
@time scaleresult = blocking(:action => scaleforgmera , :correlations => blockedcorrelations, :crystal => blockedcrystal)

scaledcrystal = scaleresult[:crystal]
scaledcorrelations = scaleresult[:correlations]
scaler = scaleresult[:transformer]

scaledcrystal = scaledcorrelations.inspace |> getcrystal

refcrystalpoints::Subset{Offset} = latticepoints(scaledcrystal)
refsamplepoints::Subset{Offset} = refcrystalpoints + c6^2 * refcrystalpoints + c6^4 * refcrystalpoints
refblockedmodes::Subset{Mode} = quantize(:pos, scaledcrystal.unitcell, 1)
refphysicalmodes::Subset{Mode} = spanoffset(refblockedmodes, refsamplepoints)

localcenter1::Offset = [1/3, 0] ∈ (refblockedmodes |> getspace)
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

emptyseedingfocks1 = (FockSpace{Region}(circularregionmodes(scaledcorrelations, pos, refphysicalmodes, 0.01)) for pos in emptyseedspos)
emptyseeds1 = sum([reduce(+ ,FockMap(v,inspace = v |> getinspace |> orderedmodes |> FockSpace,performpermute = false) for (n, (_ , v)) in enumerate(_findlocalspstates(statecorrelations = scaledcorrelations, regionfock = emptyseedingfock, 
spectrumextractpredicate = v -> true, symmetry = identitytransform(2)))) for emptyseedingfock in emptyseedingfocks1])
extendedemptyseeds1 = idmap(localfock1, localfock1)[:,emptyseeds1 |> getoutspace] * emptyseeds1
wannierizedempty1 = _localwannierization(localiso1[:empty], extendedemptyseeds1)

filledseedingfocks1 = (FockSpace{Region}(circularregionmodes(scaledcorrelations, pos, refphysicalmodes, 0.01)) for pos in filledseedspos)
filledseeds1 = sum([reduce(+ ,FockMap(v,inspace = v |> getinspace |> orderedmodes |> FockSpace,performpermute = false) for (n, (_ , v)) in enumerate(_findlocalspstates(statecorrelations = scaledcorrelations, regionfock = filledseedingfock, 
spectrumextractpredicate = v -> true, symmetry = identitytransform(2)))) for filledseedingfock in filledseedingfocks1])
extendedfilledseeds1 = idmap(localfock1, localfock1)[:,filledseeds1 |> getoutspace] * filledseeds1
wannierizedfilled1 = _localwannierization(localiso1[:filled], extendedfilledseeds1)

visualize(wannierizedempty1'*regioncorrelations(scaledcorrelations, localfock1)*wannierizedempty1|>eigspec)


emptyseedspos, filledseedspos = constructfilledandemptyseed(modebydist1[3],c6|> recenter(localcenter1))
function constructfilledandemptyseed(modewithdistlist,symmetry::AffineTransform)
    firstpos = modewithdistlist[1][2]|>getpos
    emptyseedspos = [firstpos,symmetry^2*firstpos,symmetry^4*firstpos]
    filledseedspos = [symmetry*firstpos,symmetry^3*firstpos,symmetry^5*firstpos]
    return emptyseedspos, filledseedspos
end
c6
for modewithdist in modebydist1[3]
    println(modewithdist[2]|>getpos)
end
emptyseedspos

frozenseedingfocks1 = (FockSpace{Region}(circularregionmodes(scaledcorrelations, modewithdist[2]|> getpos, refphysicalmodes, 0.01)) for modewithdist in modebydist1[3])

localcenter2::Offset = [0, 1/3] ∈ (refblockedmodes |> getspace)
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

localcenter3::Offset = [-1/3, -1/3] ∈ (refblockedmodes |> getspace)
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

# visualize((wannierizedfrozens3'*regioncorrelations(scaledcorrelations, localfock3)*wannierizedfrozens3) |> eigspec)

frozenseeds = frozenseeds1 + frozenseeds2 + frozenseeds3
wannierizedfrozens = wannierizedfrozens1 + wannierizedfrozens2 + wannierizedfrozens3

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

(frozencorrelation.inspace |> getcrystal).unitcell
visualize((frozencorrelation.inspace |> getcrystal).unitcell)
trialcrystalpoints::Subset{Offset} = latticepoints(frozencorrelation.inspace |> getcrystal)
trialsamplepoints::Subset{Offset} = trialcrystalpoints + c6^2 * trialcrystalpoints + c6^4 * trialcrystalpoints
trialblockedmodes::Subset{Mode} = quantize(:pos, (frozencorrelation.inspace |> getcrystal).unitcell, 1)
trialphysicalmodes::Subset{Mode} = spanoffset(trialblockedmodes, trialsamplepoints)
trialmodes = circularregionmodes(frozencorrelation, [1,1] ∈ (trialblockedmodes |> getspace), trialphysicalmodes, 0.3333)
trialfock = FockSpace{Region}(trialmodes)
visualize(getregion(trialfock),visualspace=euclidean(RealSpace, 2))
visualize(regioncorrelations(frozencorrelation, trialfock) |> eigspec)



modeattrs(localfock)
modeattrs(blockedcorrelations.outspace)
modeattrs(trialfock)
modeattrs(frozencorrelation.outspace)
# modeattrs(newfrozencorrelation.outspace)

# new_inspace = FockSpace(mode|>setattr(:flavor=>1) |> removeattr(:ind) |> removeattr(:dumind) for mode in frozencorrelation.inspace)
# new_crystalinspace = FockSpace(new_inspace,reflected=scaledcrystal)
# newfrozencorrelation = FockMap(frozencorrelation,inspace = new_crystalinspace,outspace = new_crystalinspace,performpermute=false)


using Distributed
procs = addprocs(5)
@info("Using $(nworkers()) threads...")

@everywhere using Plotly, SmithNormalForm, LinearAlgebra, OrderedCollections, SparseArrays, Combinatorics, DataFrames
@everywhere using Zipper

# function _crystalisometries(; localisometry::FockMap,crystalfock::FockSpace{Crystal},localregionfock::RegionFock, 
#     homemappings::Dict{Mode, Mode} = mapunitcellfock(crystalfock, localregionfock), addinspacemomentuminfo::Bool = false)

#     crystal::Crystal = getcrystal(crystalfock)
    
#     fouriermap::FockMap = fourier(blockedcorrelations|>getoutspace, localregionfock, homemappings)
#     momentumfouriers::Base.Generator = rowsubmaps(fouriermap)
#     bz::Subset{Momentum} = brillouinzone(crystal)

#     function preprocesslocalisometry(k::Momentum)::FockMap
#         if !addinspacemomentuminfo
#             return localisometry
#         end
#         inspace::FockSpace = localisometry.inspace |> orderedmodes |> setattr(:k => k) |> removeattr(:r) |> FockSpace
#         return FockMap(localisometry, inspace=inspace, performpermute=false)
#     end

#     return Dict(k => kfourier * preprocesslocalisometry(k) for (k, kfourier) in zip(bz, momentumfouriers)),fouriermap
# end

function _crystalisometries(; localisometry::FockMap, crystalfock::FockSpace{Crystal},
    addinspacemomentuminfo::Bool = false)

    crystal::Crystal = getcrystal(crystalfock)
    fouriermap::FockMap = fourier(crystalfock, localisometry|>getoutspace|>RegionFock) 
    momentumfouriers::Base.Generator = rowsubmaps(fouriermap)
    bz::Subset{Momentum} = brillouinzone(crystal)

    function preprocesslocalisometry(k::Momentum)::FockMap
        if !addinspacemomentuminfo
            return localisometry
        end
        inspace::FockSpace = localisometry.inspace |> orderedmodes |> setattr(:k => k) |> removeattr(:r) |> FockSpace
        return FockMap(localisometry, inspace=inspace, performpermute=false)
    end

    return Dict(k => kfourier * preprocesslocalisometry(k) for (k, kfourier) in zip(bz, momentumfouriers))
end

# function _crystalisometriesprime(; localisometry::FockMap, crystalfock::FockSpace{Crystal},
#     addinspacemomentuminfo::Bool = false)

#     crystal::Crystal = getcrystal(crystalfock)
#     # fouriermap::FockMap = fourier(crystalfock, localisometry.outspace|>RegionFock) / (crystal |> vol |> sqrt) seems do not need the denominator
#     fouriermap::FockMap = fourier(crystalfock, localisometry.outspace|>RegionFock) 
#     momentumfouriers::Base.Generator = rowsubmaps(fouriermap)
#     bz::Subset{Momentum} = brillouinzone(crystal)


#     function preprocesslocalisometry(k::Momentum)::FockMap
#         if !addinspacemomentuminfo
#             return localisometry
#         end
#         inspace::FockSpace = localisometry.inspace |> orderedmodes |> setattr(:offset => k) |> FockSpace
#         return FockMap(localisometry, inspace=inspace, performpermute=false)
#     end

#     return Dict(k => kfourier * preprocesslocalisometry(k) for (k, kfourier) in zip(bz, momentumfouriers)),fouriermap
# end

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
        return inmode |> setattr(:r => offset) |> setattr(:b => basis) 
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


function constructfilledandemptyseed(modewithdistlist,symmetry::AffineTransform)
    refmode = RegionFock(FockSpace(modewithdistlist[1][2]))
    emptymodes = [refmode,symmetry^2*refmode,symmetry^4*refmode]
    filledmodes = [symmetry*refmode,symmetry^3*refmode,symmetry^5*refmode]
    return emptymodes, filledmodes
end

triangular = RealSpace([sqrt(3)/2 -1/2; 0. 1.]')
kspace = convert(MomentumSpace, triangular)

pa = [1/3, 2/3] ∈ triangular
pb = [2/3, 1/3] ∈ triangular
pc = (pa + pb) / 2
spatialsnappingcalibration((pa, pb, pc))

c6 = pointgrouptransform([cos(π/3) -sin(π/3); sin(π/3) cos(π/3)])

unitcell = Subset(pa, pb)
crystal = Crystal(unitcell, [48, 48])
reciprocalhashcalibration(crystal.sizes)

modes::Subset{Mode} = quantize(unitcell, 1)|>orderedmodes
m0, m1 = members(modes)

tₙ = -1 + 0im
tₕ = 0.4im

nearestneighbor = [
    (m0, m1) => tₙ,
    (m0, m1 |> setattr(:r => [-1, 0] ∈ triangular)) => tₙ,
    (m0, m1 |> setattr(:r => [0, 1] ∈ triangular)) => tₙ]

haldane = [
    (m0, setattr(m0, :r => Point([1, 1], triangular))) => tₕ,
    (m0, setattr(m0, :r => Point([-1, 0], triangular))) => tₕ,
    (m0, setattr(m0, :r => Point([0, -1], triangular))) => tₕ,
    (m1, setattr(m1, :r => Point([1, 1], triangular))) => -tₕ,
    (m1, setattr(m1, :r => Point([-1, 0], triangular))) => -tₕ,
    (m1, setattr(m1, :r => Point([0, -1], triangular))) => -tₕ]

bonds::FockMap = bondmap([nearestneighbor..., haldane...])

energyspectrum = @time computeenergyspectrum(bonds, crystal=crystal)
energyspectrum|>visualize

groundstates::CrystalSpectrum = groundstatespectrum(energyspectrum, perunitcellfillings=1)
groundstates|>visualize

groundstateprojector = groundstates |> crystalprojector
correlations = idmap(groundstateprojector|>getoutspace) - groundstateprojector

Base.:show(io::IO, ::CrystalFockMap) = print(io, string("CrystalFockMap"))


crystalfock = correlations|>getoutspace
scale = Scale([6 0; 0 6], crystalfock|>getcrystal|>getspace)

@info("Performing blocking...")
@info("Generating blocking transformation...")
blocker = @time scale * crystalfock
@info("Performing blocking on correlations...")
blockedcorrelations = @time blocker * correlations * blocker'
blockedcrystalfock = blockedcorrelations|>getoutspace
blockedcrystal::Crystal = blockedcrystalfock|>getcrystal
blockedspace::RealSpace = blockedcrystal|>getspace

localcenter1input = [1/3, 0] ∈ blockedspace 
localcenter2input = [0, 1/3] ∈ blockedspace 
localcenter3input = [-1/3, -1/3] ∈ blockedspace 


function firststepgmera(correlations::FockMap,localcenter1,localcenter2,localcenter3,localcenter4,localcenter5)
    crystal::Crystal = getcrystal(correlations|>getinspace)
    space::RealSpace = crystal|>getspace

    function fulllocalwannierization(correlations::FockMap,localcenter::Offset,radius::Real=1/3,selectionstragedy::Function=modeselectionbycount(3))
        localregion::Region = getsphericalregion(crystal=crystal, radius=radius, metricspace=space|>orthospace)+localcenter
        localseedingfock::RegionFock = quantize(localregion,1)
        localcorr = regioncorrelations(correlations,localseedingfock)
        modebydist = groupmodesbydist(region = localregion,regionfock = localseedingfock,center = localcenter) 
        localiso = localisometries(correlations, localseedingfock, selectionstrategy=selectionstragedy)

        frozenseeds = idmap(localseedingfock, localseedingfock)[:,FockSpace(mode[2] for mode in modebydist[3])] 
        courierseeds = idmap(localseedingfock, localseedingfock)[:,FockSpace(mode[2] for mode in modebydist[1])+FockSpace(mode[2] for mode in modebydist[2])]
        wannierfrozens = _localwannierization(localiso[:filled]+localiso[:empty], frozenseeds)
        wanniercouriers = _localwannierization(localiso[:courier], courierseeds)

        return wannierfrozens, wanniercouriers,frozenseeds,courierseeds,localcorr
    end
    
    wannierizedfrozens1, wannierizedcouriers1,frozenseeds1, courierseeds1,localcorr1 = fulllocalwannierization(correlations, localcenter1)
    wannierizedfrozens2, wannierizedcouriers2,frozenseeds2, courierseeds2,localcorr2 = fulllocalwannierization(correlations, localcenter2)
    wannierizedfrozens3, wannierizedcouriers3,frozenseeds3, courierseeds3,localcorr3 = fulllocalwannierization(correlations, localcenter3)
    wannierizedfrozens4, wannierizedcouriers4,frozenseeds4, courierseeds4,localcorr4 = fulllocalwannierization(correlations, localcenter4)
    wannierizedfrozens5, wannierizedcouriers5,frozenseeds5, courierseeds5,localcorr5= fulllocalwannierization(correlations, localcenter5)

    localwannierizedfrozens = wannierizedfrozens1 + wannierizedfrozens2 + wannierizedfrozens3 + wannierizedfrozens4 + wannierizedfrozens5
    localwannierizedcouriers =  wannierizedcouriers1 + wannierizedcouriers2 + wannierizedcouriers3 + wannierizedcouriers4 + wannierizedcouriers5

    localisometry = localwannierizedfrozens+localwannierizedcouriers
    localcorr = localcorr1 + localcorr2 + localcorr3 +localcorr4 +localcorr5
    origin = [0, 0] ∈ blockedspace 
    r2unictcellfockout = FockSpace(Subset(mode for mode in localisometry |> getoutspace |> orderedmodes if (mode|>getattr(:r) == origin)))
    r2unictcellfockin = FockSpace(Subset(mode for mode in localisometry |> getoutspace |> orderedmodes if (mode|>getattr(:r) == origin)))
    unitcelllocalisometry = localisometry[r2unictcellfockin,r2unictcellfockout]
    frozenseeds = frozenseeds1 + frozenseeds2 + frozenseeds3 + frozenseeds4 + frozenseeds5
    courierseeds = courierseeds1 + courierseeds2 + courierseeds3 + courierseeds4 + courierseeds5
    localseeds = frozenseeds + courierseeds

    function globalwannierfunction(seeds, localisometry::FockMap)
        wanniercrystalisos = _crystalisometries(localisometry=localisometry, crystalfock=correlations|>getoutspace, addinspacemomentuminfo = true)

        wannierunitcell::Subset{Offset} = Subset((mode |> getattr(:b)) for mode in seeds|>getinspace |> orderedmodes)
        wanniercrystall::Crystal = Crystal(wannierunitcell, (correlations|>getoutspace|> getcrystal).sizes)
        globalwannierizedfunction::FockMap = crystaldirectsum(outcrystal = correlations|>getoutspace|>getcrystal, incrystal=wanniercrystall,wanniercrystaliso for (k,wanniercrystaliso) in wanniercrystalisos)

        return globalwannierizedfunction
    end

    # wannierfrozenisometry, fmapfrozen = globalwannierfunction(frozenseeds,localwannierizedfrozens)
    # wanniercourierisometry, fmapcourier = globalwannierfunction(courierseeds, localwannierizedcouriers)
    wannierisometry = globalwannierfunction(localseeds,unitcelllocalisometry)

    # frozencorrelations = wannierfrozenisometry' * correlations * wannierfrozenisometry
    # couriercorrelations = wanniercourierisometry' * correlations * wanniercourierisometry

    # return Dict(
    #     :frozenisometry => wannierfrozenisometry,
    #     :courierisometry => wanniercourierisometry,
    #     :frozencorrelations => frozencorrelations,
    #     :couriercorrelations => couriercorrelations,
    #     :fmapfrozen => fmapfrozen,
    #     :fmapcourier => fmapcourier)

    return Dict(
        :localwannierizedfrozens => localwannierizedfrozens,
        :localwannierizedcouriers => localwannierizedcouriers,
        :localisometry => localisometry,
        :localcorr => localcorr,
        :wannierisometry => wannierisometry)
end


localcenter1input = [1/3, 0] ∈ blockedspace 
localcenter2input = [0, 1/3] ∈ blockedspace 
localcenter3input = [1/3, 1] ∈ blockedspace 
localcenter4input = [2/3, 2/3] ∈ blockedspace 
localcenter5input = [1, 1/3] ∈ blockedspace 

# altunitvec1 = [-1/3, 1/3] ∈ blockedspace 
# altunitvec2 = [-2/3, -1/3] ∈ blockedspace 
origin = [0, 0] ∈ blockedspace 
origin == ([-0, -0] ∈ blockedspace )
firstgmerarea = firststepgmera(blockedcorrelations,localcenter1input,localcenter2input,localcenter3input,localcenter4input,localcenter5input)
visualize((firstgmerarea[:localisometry]'*firstgmerarea[:localcorr]*firstgmerarea[:localisometry])|> eigspech)
text = FockSpace(Subset(mode for mode in firstgmerarea[:localisometry] |> getoutspace |> orderedmodes if (mode|>getattr(:r) == origin)))
textp = FockSpace(Subset(mode for mode in firstgmerarea[:localisometry] |> getinspace |> orderedmodes if (mode|>getattr(:r) == origin)))

visualize(firstgmerarea[:localisometry]*firstgmerarea[:localisometry]')
firstgmerarea[:wannierizedfrozens1] 
unitcellcorr = regioncorrelations(blockedcorrelations,RegionFock(blockedcorrelations|>getoutspace|>unitcellfock))
visualize(unitcellcorr|>eigspech)
a= 1
# reflocalfock1 = RegionFock((firstgmerarea[:wannierizedcouriers1]|> getinspace) + (firstgmerarea[:wannierizedfrozens1]|> getinspace)) 
# pairmodes1 = Dict(((mode|>getattr(:b)))=>mode for mode in reflocalfock1)
# reflocalfock2 = RegionFock((firstgmerarea[:wannierizedcouriers2]|> getinspace) + (firstgmerarea[:wannierizedfrozens2]|> getinspace)) 
# pairmodes2 = Dict(((mode|>getattr(:b)))=>mode for mode in reflocalfock2)
# reflocalfock3 = RegionFock((firstgmerarea[:wannierizedcouriers3]|> getinspace) + (firstgmerarea[:wannierizedfrozens3]|> getinspace)) 
# pairmodes3 = Dict(((mode|>getattr(:b)))=>mode for mode in reflocalfock3)
# basismodes = Dict((mode|>getattr(:b)|>basispoint)=>mode for (m,mode) in blockedcorrelations|>getoutspace|>unitcellfock|>enumerate)

# reflocalfock = RegionFock(reflocalfock1 + reflocalfock2 + reflocalfock3)
# homemappings = merge(Dict(basismodes[key]=>val|>removeattr(:r) for (key,val) in pairmodes1),Dict(basismodes[key]=>val|>removeattr(:r) for (key,val) in pairmodes2),Dict(basismodes[key]=>val|>removeattr(:r) for (key,val) in pairmodes3))
# translatedmodes1wrt1 = (((mode|>getattr(:b)) + (mode|>getattr(:r))+altunitvec1)=>mode for mode in reflocalfock1)
# translatedmodes2wrt1 = (((mode|>getattr(:b)) + (mode|>getattr(:r))+altunitvec2)=>mode for mode in reflocalfock1)
# pairmodes1 = Dict(((mode|>getattr(:b)))=>mode for mode in reflocalfock1)
# pairmodes2 = (((mode|>getattr(:b)) + (mode|>getattr(:r)))=>(mode|>getattr(:b)) for mode in reflocalfock2)
# pairmodes3 = (((mode|>getattr(:b)) + (mode|>getattr(:r)))=>(mode|>getattr(:b)) for mode in reflocalfock3)

# pairmodes2wrt1 = Dict(b=> d for ((a,b),(c,d)) in zip(pairmodes2,translatedmodes1wrt1) if a==c)
# pairmodes3wrt1 = Dict(b=> d for ((a,b),(c,d)) in zip(pairmodes3,translatedmodes2wrt1) )
# basismodes = Dict((mode|>getattr(:b)|>basispoint)=>mode for (m,mode) in blockedcorrelations|>getoutspace|>unitcellfock|>enumerate)
# homemappings = merge(Dict(basismodes[key]=>val|>removeattr(:r) for (key,val) in pairmodes1),Dict(basismodes[key]=>val|>removeattr(:r) for (key,val) in pairmodes2wrt1),Dict(basismodes[key]=>val|>removeattr(:r) for (key,val) in pairmodes3wrt1))


# for (key,val) in homemappings
#     println(key|>getattr(:b),(val|>getattr(:r)),(val|>getattr(:b)))
# end

# (firstgmerarea[:wannierizedcouriers1]+firstgmerarea[:wannierizedcouriers2]+firstgmerarea[:wannierizedcouriers3])|>getoutspace|>RegionFock 
# reflocalfock
# visualize(fourier(blockedcorrelations|>getoutspace, reflocalfock, homemappings))

# wanniercrystalisos,fmap = _crystalisometries(localisometry=firstgmerarea[:localwannierizedcouriers],crystalfock=blockedcorrelations|>getoutspace,localregionfock=reflocalfock, 
#     homemappings=homemappings, addinspacemomentuminfo = true)
# wannierunitcell::Subset{Offset} = Subset((mode |> getattr(:b)) for mode in firstgmerarea[:courierseeds]|>getinspace |> orderedmodes)
# wanniercrystall::Crystal = Crystal(wannierunitcell, (correlations|>getoutspace|> getcrystal).sizes)
# globalwannierizedfunction = crystaldirectsum(outcrystal = blockedcorrelations|>getoutspace|>getcrystal, incrystal=wanniercrystall,wanniercrystaliso for (k,wanniercrystaliso) in wanniercrystalisos)
# visualize([val for (key,val) in wanniercrystalisos][51])

# wanniercrystalisosprime,fmapprime = _crystalisometriesprime(localisometry=firstgmerarea[:localwannierizedcouriers],crystalfock=blockedcorrelations|>getoutspace,localregionfock=reflocalfock, addinspacemomentuminfo = true)
# globalwannierizedfunctionprime = crystaldirectsum(outcrystal = blockedcorrelations|>getoutspace|>getcrystal, incrystal=wanniercrystall,wanniercrystaliso for (k,wanniercrystaliso) in wanniercrystalisosprime)
# visualize(fmapprime)
# visualize(fmap)
firstgmerarea[:wannierisometry]'*firstgmerarea[:wannierisometry]

refcorrelations = firstgmerarea[:wannierisometry]'*blockedcorrelations*firstgmerarea[:wannierisometry]
visualize((firstgmerarea[:wannierisometry]'*firstgmerarea[:wannierisometry])|>crystalspectrum)
visualize((refcorrelations|>getinspace |> getcrystal).unitcell)
trialcrystal::Crystal = refcorrelations|> getoutspace|>getcrystal
trialspace::RealSpace = trialcrystal|>getspace

triallocalregion::Region = getsphericalregion(crystal=trialcrystal, radius=1/3, metricspace=trialspace|>orthospace)
visualize(triallocalregion)
triallocalseedingfock::RegionFock = quantize(triallocalregion, 1)
visualize(regioncorrelations(refcorrelations, triallocalseedingfock)|>eigspech)

visualize(regioncorrelations(blockedcorrelations, RegionFock(blockedcorrelations|>getoutspace|>unitcellfock))|>eigspech)

localcenter1input = [1/3, 0] ∈ blockedspace 
localcenter2input = [1/3, 1/3] ∈ blockedspace 
localcenter3input = [0, 0] ∈ blockedspace 
localregion1::Region = getsphericalregion(crystal=blockedcrystal, radius=1/3, metricspace=blockedspace|>orthospace)+localcenter1input
localregion2::Region = getsphericalregion(crystal=blockedcrystal, radius=1/3, metricspace=blockedspace|>orthospace)+localcenter2input
localregion3::Region = getsphericalregion(crystal=blockedcrystal, radius=1/3, metricspace=blockedspace|>orthospace)+localcenter3input
localseedingfock1::RegionFock = quantize(localregion1,1)
localseedingfock2::RegionFock = quantize(localregion2,1)
localseedingfock3::RegionFock = quantize(localregion3,1)
visualize(regioncorrelations(blockedcorrelations, localseedingfock2)|>eigspech)

testregion::Region = getsphericalregion(crystal=blockedcrystal, radius=1/3, metricspace=blockedspace|>orthospace)+ localcenter1input
testseedingfock::RegionFock = quantize(testregion,1)
visualize(regioncorrelations(blockedcorrelations,testseedingfock)|>eigspech)

visualize(regioncorrelations(blockedcorrelations, RegionFock(localseedingfock1+localseedingfock2+localseedingfock3))|>eigspech)
refcorrelations|>getinspace|>unitcellfock|>modeattrs
triallocalseedingfock|>unitcellfock|>modeattrs
rulemapunitcellfock 



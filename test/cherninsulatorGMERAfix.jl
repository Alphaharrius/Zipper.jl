using Distributed
procs = addprocs(5)
@info("Using $(nworkers()) threads...")

@everywhere using Plotly, SmithNormalForm, LinearAlgebra, OrderedCollections, SparseArrays, Combinatorics, DataFrames
@everywhere using Zipper

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

# function _crystalisometries(; localisometry::FockMap, crystalfock::FockSpace{Crystal},
#     addinspacemomentuminfo::Bool = false)::Dict{Momentum, FockMap}

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

#     return Dict(k => kfourier * preprocesslocalisometry(k) for (k, kfourier) in zip(bz, momentumfouriers))
# end

# function _crystalisometries(; localisometry::FockMap, crystalfock::FockSpace{Crystal},
#     addinspacemomentuminfo::Bool = false)::Dict{Momentum, FockMap}

#     crystal::Crystal = getcrystal(crystalfock)
#     fouriermap::FockMap = fourier(crystalfock, localisometry.outspace|>RegionFock) / (crystal |> vol |> sqrt)
#     momentumfouriers::Base.Generator = rowsubmaps(fouriermap)
#     bz::Subset{Momentum} = brillouinzone(crystal)

#     function preprocesslocalisometry(k::Momentum)::FockMap
#         if !addinspacemomentuminfo
#             return localisometry
#         end
#         inspace::FockSpace = localisometry|>getinspace|>orderedmodes|>setattr(:k=>k)|>FockSpace
#         return FockMap(localisometry, inspace=inspace, performpermute=false)
#     end

#     return Dict(k => kfourier * preprocesslocalisometry(k) for (k, kfourier) in zip(bz, momentumfouriers))
# end

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


localregion::Region = getsphericalregion(crystal=blockedcrystal, radius=1/3, metricspace=blockedspace|>orthospace)+localcenter
localseedingfock::RegionFock = quantize(localregion, 1)
modebydist = groupmodesbydist(region = localregion,regionfock = localseedingfock,center = localcenter) 
modebydist[3]
# emptymodes, filledmodes = constructfilledandemptyseed(modebydist[3],c6|> recenter(localcenter))
emptymodes

localiso = localisometries(blockedcorrelations,localseedingfock,selectionstrategy=modeselectionbycount(3))
frozenseeds = idmap(localseedingfock, localseedingfock)[:,FockSpace(mode[2] for mode in modebydist[3])] 
courierseeds = idmap(localseedingfock, localseedingfock)[:,FockSpace(mode[2] for mode in modebydist[1])+FockSpace(mode[2] for mode in modebydist[2])] 
_localwannierization(localiso[:filled]+localiso[:empty], frozenseeds)
_localwannierization(localiso[:courier], courierseeds)

idmap(localseedingfock, localseedingfock)*(localiso[:filled]+localiso[:empty])


function firststepgmera(correlations::FockMap,localcenter1,localcenter2,localcenter3)
    crystal::Crystal = getcrystal(correlations|>getinspace)
    space::RealSpace = crystal|>getspace

    function fulllocalwannierization(correlations::FockMap,localcenter::Offset,radius::Real=1/3,selectionstragedy::Function=modeselectionbycount(3))
        localregion::Region = getsphericalregion(crystal=crystal, radius=radius, metricspace=space|>orthospace)+localcenter
        localseedingfock::RegionFock = quantize(localregion,1)
        modebydist = groupmodesbydist(region = localregion,regionfock = localseedingfock,center = localcenter) 
        localiso = localisometries(correlations, localseedingfock, selectionstrategy=selectionstragedy)

        frozenseeds = idmap(localseedingfock, localseedingfock)[:,FockSpace(mode[2] for mode in modebydist[3])] 
        courierseeds = idmap(localseedingfock, localseedingfock)[:,FockSpace(mode[2] for mode in modebydist[1])+FockSpace(mode[2] for mode in modebydist[2])]
        wannierfrozens = _localwannierization(localiso[:filled]+localiso[:empty], frozenseeds)
        wanniercouriers = _localwannierization(localiso[:courier], courierseeds)

        return wannierfrozens, wanniercouriers,frozenseeds,courierseeds
    end
    
    wannierizedfrozens1, wannierizedcouriers1,frozenseeds1, courierseeds1 = fulllocalwannierization(correlations, localcenter1)
    wannierizedfrozens2, wannierizedcouriers2,frozenseeds2, courierseeds2 = fulllocalwannierization(correlations, localcenter2)
    wannierizedfrozens3, wannierizedcouriers3,frozenseeds3, courierseeds3 = fulllocalwannierization(correlations, localcenter3)

    localwannierizedfrozens = wannierizedfrozens1 + wannierizedfrozens2 + wannierizedfrozens3
    localwannierizedcouriers =  wannierizedcouriers1 + wannierizedcouriers2 + wannierizedcouriers3

    frozenseeds = frozenseeds1 + frozenseeds2 + frozenseeds3
    courierseeds = courierseeds1 + courierseeds2 + courierseeds3

    function globalwannierfunction(seeds, localisometry::FockMap)
        wanniercrystalisos = crystalisometries(localisometry=localisometry, crystalfock=correlations|>getoutspace, addinspacemomentuminfo = true)
        wannierunitcell::Subset{Offset} = Subset((mode |> getattr(:b)) for mode in seeds|>getinspace |> orderedmodes)
        wanniercrystall::Crystal = Crystal(wannierunitcell, (correlations|>getoutspace|> getcrystal).sizes)
        globalwannierizedfunction::FockMap = crystaldirectsum(outcrystal = correlations|>getoutspace|>getcrystal, incrystal=wanniercrystall,wanniercrystaliso for (k,wanniercrystaliso) in wanniercrystalisos)

        return globalwannierizedfunction
    end

    wannierfrozenisometry = globalwannierfunction(frozenseeds,localwannierizedfrozens)
    wanniercourierisometry = globalwannierfunction(courierseeds, localwannierizedcouriers)

    frozencorrelations = wannierfrozenisometry' * correlations * wannierfrozenisometry
    couriercorrelations = wanniercourierisometry' * correlations * wanniercourierisometry

    return Dict(
        :frozenisometry => wannierfrozenisometry,
        :courierisometry => wanniercourierisometry,
        :frozencorrelations => frozencorrelations,
        :couriercorrelations => couriercorrelations)
end


localcenter1input = [1/3, 0] ∈ blockedspace 
localcenter2input = [0, 1/3] ∈ blockedspace 
localcenter3input = [-1/3, -1/3] ∈ blockedspace 

firstgmerarea = firststepgmera(blockedcorrelations,localcenter1input,localcenter2input,localcenter3input)

refcorrelations = firstgmerarea[:couriercorrelations]
visualize((refcorrelations|>getinspace |> getcrystal).unitcell)

Subset(mode|>getpos for mode in refcorrelations|>getoutspace |>orderedmodes)

trialcrystal::Crystal = refcorrelations|>getoutspace|>getcrystal
refcorrelations|>getoutspace |> orderedmodes

trialspace::RealSpace = trialcrystal|>getspace

triallocalregion::Region = getsphericalregion(crystal=trialcrystal, radius=1/3, metricspace=trialspace|>orthospace)
visualize(triallocalregion)


triallocalseedingfock::RegionFock = quantize(triallocalregion, 1)
refcorrelations|>getoutspace|>modeattrs
triallocalseedingfock|>modeattrs
regioncorrelations(refcorrelations, triallocalseedingfock)





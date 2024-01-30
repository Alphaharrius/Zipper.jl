using Distributed
procs = addprocs(5)
@info("Using $(nworkers()) threads...")

@everywhere using Plotly, SmithNormalForm, LinearAlgebra, OrderedCollections, SparseArrays, Combinatorics, DataFrames
@everywhere using Zipper

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
localcenter3input = [1/3, 1] ∈ blockedspace 
localcenter4input = [2/3, 2/3] ∈ blockedspace 
localcenter5input = [1, 1/3] ∈ blockedspace 

function fulllocalwannierization(correlations::FockMap,localcenter::Offset,radius::Real=1/3,selectionstragedy::Function=modeselectionbycount(3))
    localregion::Region = getsphericalregion(crystal=blockedcrystal, radius=radius, metricspace=blockedspace|>orthospace)+localcenter
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

function globalwannierfunction(seeds, localisometry::FockMap)
    wanniercrystalisos = _crystalisometries(localisometry=localisometry, crystalfock=correlations|>getoutspace, addinspacemomentuminfo = true)

    wannierunitcell::Subset{Offset} = Subset((mode |> getattr(:b)) for mode in seeds|>getinspace |> orderedmodes)
    wanniercrystall::Crystal = Crystal(wannierunitcell, (correlations|>getoutspace|> getcrystal).sizes)
    globalwannierizedfunction::FockMap = crystaldirectsum(outcrystal = correlations|>getoutspace|>getcrystal, incrystal=wanniercrystall,wanniercrystaliso for (k,wanniercrystaliso) in wanniercrystalisos)

    return globalwannierizedfunction
end

wannierizedfrozens1, wannierizedcouriers1,frozenseeds1, courierseeds1,localcorr1 = fulllocalwannierization(blockedcorrelations, localcenter1input)
wannierizedfrozens2, wannierizedcouriers2,frozenseeds2, courierseeds2,localcorr2 = fulllocalwannierization(blockedcorrelations, localcenter2input)
wannierizedfrozens3, wannierizedcouriers3,frozenseeds3, courierseeds3,localcorr3 = fulllocalwannierization(blockedcorrelations, localcenter3input)
wannierizedfrozens4, wannierizedcouriers4,frozenseeds4, courierseeds4,localcorr4 = fulllocalwannierization(blockedcorrelations, localcenter4input)
wannierizedfrozens5, wannierizedcouriers5,frozenseeds5, courierseeds5,localcorr5 = fulllocalwannierization(blockedcorrelations, localcenter5input)

extendedwannierizedfrozens = wannierizedfrozens1 + wannierizedfrozens2 + wannierizedfrozens3 + wannierizedfrozens4 + wannierizedfrozens5
extendedwannierizedcouriers =  wannierizedcouriers1 + wannierizedcouriers2 + wannierizedcouriers3 + wannierizedcouriers4 + wannierizedcouriers5
extendedwannierlocaliso = extendedwannierizedfrozens + extendedwannierizedcouriers
extendedcourierseeds = courierseeds1 + courierseeds2 + courierseeds3 + courierseeds4 + courierseeds5
extendedfrozenseeds = frozenseeds1 + frozenseeds2 + frozenseeds3 + frozenseeds4 + frozenseeds5
extendedlocalseeds = extendedcourierseeds + extendedfrozenseeds

# origin = [0, 0] ∈ blockedspace 
# r2unictcellfock = FockSpace(Subset(mode for mode in extendedlocalseeds |> getoutspace |> orderedmodes if (mode|>getattr(:r) == origin)))

# extendedlocalcorr = regioncorrelations(blockedcorrelations,RegionFock((extendedwannierlocaliso)|>getoutspace))
# rotextendedlocalcorr = extendedwannierlocaliso[:,r2unictcellfock]'*extendedlocalcorr*extendedwannierlocaliso[:,r2unictcellfock]
# visualize(rotextendedlocalcorr|>eigspech)
# extendedwannierlocaliso[:,r2unictcellfock] |> getoutspace
# wanniercrystalisos = _crystalisometries(localisometry=extendedwannierlocaliso[:,r2unictcellfock], crystalfock=blockedcorrelations|>getoutspace, addinspacemomentuminfo = true)
# wannierunitcell::Subset{Offset} = Subset((mode |> getattr(:b)) for mode in extendedlocalseeds|>getinspace |> orderedmodes)
# wanniercrystall::Crystal = Crystal(wannierunitcell, (blockedcorrelations|>getoutspace|> getcrystal).sizes)
# globalwannierizedfunction::FockMap = crystaldirectsum(outcrystal = blockedcorrelations|>getoutspace|>getcrystal, incrystal=wanniercrystall,wanniercrystaliso for (k,wanniercrystaliso) in wanniercrystalisos)


origin = [0, 0] ∈ blockedspace 
r2unictcellfock = FockSpace(Subset(mode for mode in extendedcourierseeds |> getinspace |> orderedmodes if (mode|>getattr(:r) == origin)))

extendedlocalcorr = regioncorrelations(blockedcorrelations,RegionFock((extendedwannierizedcouriers)|>getoutspace))
extendedwannierizedcouriers[:,r2unictcellfock]
rotextendedlocalcorr = extendedwannierizedcouriers[:,r2unictcellfock]'*extendedlocalcorr*extendedwannierizedcouriers[:,r2unictcellfock]
visualize(rotextendedlocalcorr|>eigspech)
# extendedwannierlocaliso[:,r2unictcellfock] |> getoutspace
wanniercrystalisos = _crystalisometries(localisometry=extendedwannierizedcouriers[:,r2unictcellfock], crystalfock=blockedcorrelations|>getoutspace, addinspacemomentuminfo = true)
wannierunitcell::Subset{Offset} = Subset((mode |> getattr(:b)) for mode in extendedcourierseeds|>getinspace |> orderedmodes)
wanniercrystall::Crystal = Crystal(wannierunitcell, (blockedcorrelations|>getoutspace|> getcrystal).sizes)
globalwannierizedfunction::FockMap = crystaldirectsum(outcrystal = blockedcorrelations|>getoutspace|>getcrystal, incrystal=wanniercrystall,wanniercrystaliso for (k,wanniercrystaliso) in wanniercrystalisos)

refcorrelations = globalwannierizedfunction'*blockedcorrelations*globalwannierizedfunction
visualize(refcorrelations|>crystalspectrum)
visualize((refcorrelations|>getinspace |> getcrystal).unitcell)
trialcrystal::Crystal = refcorrelations|> getoutspace|>getcrystal
trialspace::RealSpace = trialcrystal|>getspace

shift = [2/3, 0] ∈ blockedspace 
triallocalregion::Region = getsphericalregionwifcenter(center = origin,crystal=trialcrystal, radius=1/3, metricspace=trialspace|>orthospace) 
visualize(triallocalregion)
triallocalseedingfock::RegionFock = quantize(triallocalregion, 1)
visualize(regioncorrelations(refcorrelations, triallocalseedingfock)|>eigspech)


localregion::Region = getsphericalregionwifcenter(center = localcenter2input,crystal=blockedcrystal, radius=1/3, metricspace=blockedspace|>orthospace)
visualize(localregion)
localseedingfock::RegionFock = quantize(localregion,1)
visualize(regioncorrelations(blockedcorrelations,localseedingfock )|>eigspech)

function getsphericalregionwifcenter(;center, crystal::Crystal, radius::Real, metricspace::RealSpace)
    generatingradius::Integer = ceil(Int, radius * 1.5) # Multiply by 1.5 to ensure all unitcell points fits.
    generatinglength::Integer = generatingradius * 2
    generatingcrystal::Crystal = Crystal(crystal|>getunitcell, [generatinglength, generatinglength])
    crystalregion::Region = generatingcrystal|>sitepoints
    centeredregion::Region = crystalregion - (crystalregion|>getcenter) 
    return Subset(point for point in centeredregion if norm(metricspace * (point-center)) <= radius)
end
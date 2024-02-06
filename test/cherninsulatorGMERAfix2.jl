using Distributed
procs = addprocs(5)
@info("Using $(nworkers()) threads...")

@everywhere using Plotly, SmithNormalForm, LinearAlgebra, OrderedCollections, SparseArrays, Combinatorics, DataFrames
@everywhere using Zipper

function _crystalisometries(; localisometry::FockMap, crystalfock::CrystalFock,
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

# Maybe wrong if the offset is too far away
function getsphericalregionwifcenter(;center, crystal::Crystal, radius::Real, metricspace::RealSpace)
    generatingradius::Integer = ceil(Int, radius * 1.5) # Multiply by 1.5 to ensure all unitcell points fits.
    generatinglength::Integer = generatingradius * 2
    generatingcrystal::Crystal = Crystal(crystal|>getunitcell, [generatinglength,generatinglength])
    crystalregion::Region = generatingcrystal|>sitepoints
    centeredregion::Region = crystalregion - (crystalregion|>getcenter) 
    return Subset(point for point in crystalregion if norm(metricspace * (point-center)) <= radius)
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

    spatialinspace::RegionFock = RegionFock( _spatialinmode(fockmap[:, m],i) for (i,m) in fockmap |> getinspace |> enumerate)
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
    posmodedict = Dict(modepair[2]|>getpos => modepair[2] for modepair in modewithdistlist)
    refmode = modewithdistlist[1][2]
    c6poslist = [(symmetry^r)*(refmode|>getpos) for r in range(0,5)]
    emptyposlist = [c6poslist[1],c6poslist[3],c6poslist[5]]
    filledposlist = [c6poslist[2],c6poslist[4],c6poslist[6]]
    emptymodes = [posmodedict[emptypos] for emptypos in emptyposlist]
    filledmodes = [posmodedict[filledpos] for filledpos in filledposlist]
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

H = energyspectrum |> FockMap

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

H
blockingmap = FockMap(blocker)
courierisometry = FockMap(firstgmerarea[:courierisometry])
frozenisometry = FockMap(firstgmerarea[:frozenisometry])
visualize((courierisometry'*blockingmap*H*blockingmap'*courierisometry)|>crystalspectrum)
visualize((frozenisometry'*blockingmap*H*blockingmap'*frozenisometry)|>crystalspectrum)
visualize((blockingmap*H*blockingmap')|>crystalspectrum)
visualize((courierisometry'*blockingmap*H*blockingmap'*courierisometry+frozenisometry'*blockingmap*H*blockingmap'*frozenisometry)|>crystalspectrum)
courierisometry'*blockingmap*H*blockingmap'*courierisometry+frozenisometry'*blockingmap*H*blockingmap'*frozenisometry
frozenisometry'*blockingmap*H*blockingmap'*frozenisometry


localcenter1input = [1/3, 0] ∈ blockedspace 
localcenter2input = [0, 1/3] ∈ blockedspace 
localcenter3input = [1/3, 1] ∈ blockedspace 
localcenter4input = [2/3, 2/3] ∈ blockedspace 
localcenter5input = [1, 1/3] ∈ blockedspace 

function firststepgmera(correlations::FockMap,localcenter1,localcenter2,localcenter3,localcenter4,localcenter5)
    crystal::Crystal = getcrystal(correlations|>getinspace)
    space::RealSpace = crystal|>getspace

    function fulllocalwannierization(correlations::FockMap,localcenter::Offset,radius::Real=1/3,selectionstragedy::Function=modeselectionbycount(3))
        localregion::Region = getsphericalregionwifcenter(center=localcenter,crystal=crystal, radius=radius, metricspace=space|>orthospace)
        localseedingfock::RegionFock = quantize(localregion,1)
        modebydist = groupmodesbydist(region = localregion,regionfock = localseedingfock,center = localcenter) 
        (emptyseedslist,filledseedslist) = constructfilledandemptyseed(modebydist[3],c6|>recenter(localcenter))
        localiso = localisometries(correlations, localseedingfock, selectionstrategy=selectionstragedy)

        frozenseeds = idmap(localseedingfock, localseedingfock)[:,FockSpace(mode[2] for mode in modebydist[3])] 
        emptyseeds = idmap(localseedingfock, localseedingfock)[:,FockSpace(mode for mode in emptyseedslist)] 
        filledseeds = idmap(localseedingfock, localseedingfock)[:,FockSpace(mode for mode in filledseedslist)] 

        courierseeds = idmap(localseedingfock, localseedingfock)[:,FockSpace(mode[2] for mode in modebydist[1])+FockSpace(mode[2] for mode in modebydist[2])]
        wannierfrozens = _localwannierization(localiso[:filled]+localiso[:empty], frozenseeds)
        wannieremptys = _localwannierization(localiso[:empty], emptyseeds)
        wannierfilleds = _localwannierization(localiso[:filled], filledseeds)
        wanniercouriers = _localwannierization(localiso[:courier], courierseeds)

        return wannierfrozens,wannieremptys,wannierfilleds, wanniercouriers,frozenseeds,emptyseeds,filledseeds,courierseeds
    end
    
    wannierizedfrozens1,wannierizedemptys1,wannierizedfilleds1, wannierizedcouriers1,frozenseeds1,emptyseeds1,filledseeds1,courierseeds1 = fulllocalwannierization(correlations, localcenter1)
    wannierizedfrozens2,wannierizedemptys2,wannierizedfilleds2, wannierizedcouriers2,frozenseeds2,emptyseeds2,filledseeds2,courierseeds2 = fulllocalwannierization(correlations, localcenter2)
    wannierizedfrozens3,wannierizedemptys3,wannierizedfilleds3, wannierizedcouriers3,frozenseeds3,emptyseeds3,filledseeds3,courierseeds3 = fulllocalwannierization(correlations, localcenter3)
    wannierizedfrozens4,wannierizedemptys4,wannierizedfilleds4, wannierizedcouriers4,frozenseeds4,emptyseeds4,filledseeds4,courierseeds4 = fulllocalwannierization(correlations, localcenter4)
    wannierizedfrozens5,wannierizedemptys5,wannierizedfilleds5, wannierizedcouriers5,frozenseeds5,emptyseeds5,filledseeds5,courierseeds5 = fulllocalwannierization(correlations, localcenter5)

    extendedwannierizedfrozens = wannierizedfrozens1 + wannierizedfrozens2 + wannierizedfrozens3 + wannierizedfrozens4 + wannierizedfrozens5
    extendedwannierizedemptys =  wannierizedemptys1 + wannierizedemptys2 + wannierizedemptys3 + wannierizedemptys4 + wannierizedemptys5
    extendedwannierizedfilleds =  wannierizedfilleds1 + wannierizedfilleds2 + wannierizedfilleds3 + wannierizedfilleds4 + wannierizedfilleds5
    extendedwannierizedcouriers =  wannierizedcouriers1 + wannierizedcouriers2 + wannierizedcouriers3 + wannierizedcouriers4 + wannierizedcouriers5

    extendedcourierseeds = courierseeds1 + courierseeds2 + courierseeds3 + courierseeds4 + courierseeds5
    extendedemptyseeds = emptyseeds1 + emptyseeds2 + emptyseeds3 + emptyseeds4 + emptyseeds5
    extendedfilledseeds = filledseeds1 + filledseeds2 + filledseeds3 + filledseeds4 + filledseeds5
    extendedfrozenseeds = frozenseeds1 + frozenseeds2 + frozenseeds3 + frozenseeds4 + frozenseeds5

    origin = [0, 0] ∈ blockedspace 

    resunictcellfockcourier = FockSpace(Subset(mode for mode in extendedcourierseeds |> getinspace |> orderedmodes if (mode|>getattr(:r) == origin)))
    resunictcellfockfrozen = FockSpace(Subset(mode for mode in extendedfrozenseeds |> getinspace |> orderedmodes if (mode|>getattr(:r) == origin)))
    resunictcellfockempty = FockSpace(Subset(mode for mode in extendedemptyseeds |> getinspace |> orderedmodes if (mode|>getattr(:r) == origin)))
    resunictcellfockfilled = FockSpace(Subset(mode for mode in extendedfilledseeds |> getinspace |> orderedmodes if (mode|>getattr(:r) == origin)))

    function globalwannierfunction(seeds, localisometry::FockMap)
        wanniercrystalisos = _crystalisometries(localisometry=localisometry, crystalfock=correlations|>getoutspace, addinspacemomentuminfo = true)

        wannierunitcell::Subset{Offset} = Subset((mode |> getattr(:b)) for mode in seeds|>getinspace |> orderedmodes)
        wanniercrystall::Crystal = Crystal(wannierunitcell, (correlations|>getoutspace|> getcrystal).sizes)
        globalwannierizedfunction::FockMap = crystaldirectsum(outcrystal = correlations|>getoutspace|>getcrystal, incrystal=wanniercrystall,wanniercrystaliso for (k,wanniercrystaliso) in wanniercrystalisos)

        return globalwannierizedfunction
    end

    wanniercourierisometry = globalwannierfunction(extendedcourierseeds,extendedwannierizedcouriers[:,resunictcellfockcourier])
    wannierfrozenisometry = globalwannierfunction(extendedfrozenseeds,extendedwannierizedfrozens[:,resunictcellfockfrozen])
    wannieremptyisometry = globalwannierfunction(extendedemptyseeds,extendedwannierizedemptys[:,resunictcellfockempty])
    wannierfilledisometry = globalwannierfunction(extendedfilledseeds,extendedwannierizedfilleds[:,resunictcellfockfilled])

    couriercorrelations = wanniercourierisometry' * correlations * wanniercourierisometry
    couriercorrelationspectrum = couriercorrelations |> crystalspectrum
    purifiedcorrelationspectrum = couriercorrelationspectrum |> roundingpurification
    purifiedcouriercorrelations = purifiedcorrelationspectrum |> CrystalFockMap

    frozencorrelations = wannierfrozenisometry' * correlations * wannierfrozenisometry
    frozencorrelationspectrum = frozencorrelations |> crystalspectrum
    purifiedcorrelationspectrum = frozencorrelationspectrum |> roundingpurification
    purifiedfrozencorrelations = purifiedcorrelationspectrum |> CrystalFockMap

    emptycorrelations = wannieremptyisometry' * correlations * wannieremptyisometry
    emptycorrelationspectrum = emptycorrelations |> crystalspectrum
    purifiedcorrelationspectrum = emptycorrelationspectrum |> roundingpurification
    purifiedemptycorrelations = purifiedcorrelationspectrum |> CrystalFockMap

    filledcorrelations = wannierfilledisometry' * correlations * wannierfilledisometry
    filledcorrelationspectrum = filledcorrelations |> crystalspectrum
    purifiedcorrelationspectrum = filledcorrelationspectrum |> roundingpurification
    purifiedfilledcorrelations = purifiedcorrelationspectrum |> CrystalFockMap

    return Dict(
        :frozenisometry => wannierfrozenisometry,
        :emptyisometry => wannieremptyisometry,
        :filledisometry => wannierfilledisometry,
        :courierisometry => wanniercourierisometry,
        :frozencorrelations => frozencorrelations,
        :emptycorrelations => emptycorrelations,
        :filledcorrelations => filledcorrelations,
        :couriercorrelations => couriercorrelations,
        :purifiedfrozencorrelations => purifiedfrozencorrelations,
        :purifiedfilledcorrelations => purifiedfilledcorrelations,
        :purifiedemptycorrelations => purifiedemptycorrelations,
        :purifiedcouriercorrelations => purifiedcouriercorrelations)
end

firstgmerarea = firststepgmera(blockedcorrelations,localcenter1input,localcenter2input,localcenter3input,localcenter4input,localcenter5input)

testregion = getsphericalregionwifcenter(center=localcenter3input,crystal=blockedcrystal, radius=1/3, metricspace=blockedspace|>orthospace)
localseedingfock::RegionFock = quantize(localregion,1)

refcorrelations = firstgmerarea[:purifiedemptycorrelations]
visualize(refcorrelations|>crystalspectrum)
visualize((refcorrelations|>getinspace |> getcrystal).unitcell)
trialcrystal::Crystal = refcorrelations|> getoutspace|>getcrystal
trialspace::RealSpace = trialcrystal|>getspace

shift = [2/3, 0] ∈ blockedspace 
triallocalregion::Region,fullregion = getsphericalregionwifcenter(center = localcenter3input,crystal=trialcrystal, radius=1/3, metricspace=trialspace|>orthospace) 
visualize(triallocalregion)
visualize(fullregion)
triallocalseedingfock::RegionFock = quantize(triallocalregion, 1)
visualize(regioncorrelations(refcorrelations, triallocalseedingfock)|>eigspech)


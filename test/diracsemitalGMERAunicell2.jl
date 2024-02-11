using Plotly, SmithNormalForm, LinearAlgebra, OrderedCollections, SparseArrays, Combinatorics, DataFrames
using Zipper

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
    samedistancethreshold::Int = 3)

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
    generatingcrystal::Crystal = Crystal(crystal|>getunitcell, 4*(crystal|>size))
    crystalregion::Region = generatingcrystal|>sitepoints
    centeredregion::Region = crystalregion - (crystalregion|>getcenter) 
    return Subset(point for point in centeredregion if norm(metricspace * (point-center)) <= radius)
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
c3 = c6^2

unitcell = Subset(pa, pb)
crystal = Crystal(unitcell, [24, 24])
reciprocalhashcalibration(crystal.sizes)

modes::Subset{Mode} = quantize(unitcell, 1)|>orderedmodes
m0, m1 = members(modes)

tₙ = ComplexF64(-1.)
bonds::FockMap = bondmap([
    (m0, m1) => tₙ,
    (m0, m1 |> setattr(:r => [-1, 0] ∈ triangular)) => tₙ,
    (m0, m1 |> setattr(:r => [0, 1] ∈ triangular)) => tₙ])

energyspectrum = computeenergyspectrum(bonds, crystal=crystal)
energyspectrum |> visualize

H = energyspectrum |> FockMap

groundstates::CrystalSpectrum = groundstatespectrum(energyspectrum, perunitcellfillings=1)
groundstates |> visualize
groundstateprojector = groundstates |> crystalprojector

C = idmap(groundstateprojector|>getoutspace) - groundstateprojector

correlations = C
crystalfock = correlations|>getoutspace

localcenter = [1/2,0] ∈ blockedspace
scale = Scale([4 0; 0 4], crystalfock|>getcrystal|>getspace)
@info("Performing blocking...")
@info("Generating blocking transformation...")
blocker = @time scale * crystalfock
@info("Performing blocking on correlations...")
blockedcorrelations = @time blocker * correlations * blocker'
blockedcrystalfock = blockedcorrelations|>getoutspace
blockedcrystal::Crystal = blockedcrystalfock|>getcrystal
blockedspace::RealSpace = blockedcrystal|>getspace

shift = [1/2,1/2] ∈ blockedspace
localregion::Region =  Subset(pt+localcenter for pt in (blockedcrystal|>getunitcell))
localseedingfock::RegionFock = quantize(localregion,1)
modebydist = groupmodesbydist(region = localregion,regionfock = localseedingfock,center = localcenter+shift) 
(emptyseedslist,filledseedslist) = constructfilledandemptyseed(modebydist[6],c6|>recenter(localcenter+shift))
localiso = localisometries(blockedcorrelations, localseedingfock, selectionstrategy=modeselectionbycount(1))
visualize(regioncorrelations(blockedcorrelations, localseedingfock)|>eigspec)
emptyseeds = idmap(localseedingfock, localseedingfock)[:,FockSpace(emptyseedslist[1])] 
filledseeds = idmap(localseedingfock, localseedingfock)[:,FockSpace(filledseedslist[1])] 
courierseeds = idmap(localseedingfock, localseedingfock)[:,sum([FockSpace(mode[2] for mode in modebydist[i]) for i in range(1,5)])]

wannieremptys = _localwannierization(localiso[:empty], emptyseeds)
wannierfilleds = _localwannierization(localiso[:filled], filledseeds)
wanniercouriers = _localwannierization(localiso[:courier], courierseeds)

function gmera(correlations::FockMap, step::Int)
    crystalfock = correlations|>getoutspace
    if step == 1
        scale = Scale([4 0; 0 4], crystalfock|>getcrystal|>getspace)
        
    else
        scale = Scale([2 0; 0 2], crystalfock|>getcrystal|>getspace)
    end
    @info("Performing blocking...")
    @info("Generating blocking transformation...")
    blocker = @time scale * crystalfock
    @info("Performing blocking on correlations...")
    blockedcorrelations = @time blocker * correlations * blocker'
    blockedcrystalfock = blockedcorrelations|>getoutspace
    blockedcrystal::Crystal = blockedcrystalfock|>getcrystal
    blockedspace::RealSpace = blockedcrystal|>getspace

    function gmera1ststep(correlations::FockMap,localcenter1,localcenter2)
        crystal::Crystal = getcrystal(correlations|>getinspace)

        function fulllocalwannierization(correlations::FockMap,localcenter::Offset,selectionstragedy::Function=modeselectionbycount(3))
            shift = [1/2,1/2] ∈ blockedspace
            localregion::Region =  Subset(pt+localcenter for pt in (crystal|>getunitcell))
            localseedingfock::RegionFock = quantize(localregion,1)
            modebydist = groupmodesbydist(region = localregion,regionfock = localseedingfock,center = localcenter+shift) 
            (emptyseedslist,filledseedslist) = constructfilledandemptyseed(modebydist[6],c6|>recenter(localcenter+shift))
            localiso = localisometries(correlations, localseedingfock, selectionstrategy=selectionstragedy)
    
            emptyseeds = idmap(localseedingfock, localseedingfock)[:,FockSpace(mode for mode in emptyseedslist)] 
            filledseeds = idmap(localseedingfock, localseedingfock)[:,FockSpace(mode for mode in filledseedslist)] 
            courierseeds = idmap(localseedingfock, localseedingfock)[:,sum([FockSpace(mode[2] for mode in modebydist[i]) for i in range(1,5)])]

            wannieremptys = _localwannierization(localiso[:empty], emptyseeds)
            wannierfilleds = _localwannierization(localiso[:filled], filledseeds)
            wanniercouriers = _localwannierization(localiso[:courier], courierseeds)
    
            return wannieremptys,wannierfilleds, wanniercouriers,emptyseeds,filledseeds,courierseeds
        end
        
        wannierizedemptys1,wannierizedfilleds1, wannierizedcouriers1,emptyseeds1,filledseeds1,courierseeds1 = fulllocalwannierization(correlations, localcenter1)
        wannierizedemptys2,wannierizedfilleds2, wannierizedcouriers2,emptyseeds2,filledseeds2,courierseeds2 = fulllocalwannierization(correlations, localcenter2)
    
        extendedwannierizedemptys =  wannierizedemptys1 + wannierizedemptys2 
        extendedwannierizedfilleds =  wannierizedfilleds1 + wannierizedfilleds2
        extendedwannierizedcouriers =  wannierizedcouriers1 + wannierizedcouriers2
    
        extendedcourierseeds = courierseeds1 + courierseeds2 
        extendedemptyseeds = emptyseeds1 + emptyseeds2
        extendedfilledseeds = filledseeds1 + filledseeds2 
        
        origin = [0, 0] ∈ blockedspace 
        resunictcellfockcourier = FockSpace(Subset(mode for mode in extendedcourierseeds |> getinspace |> orderedmodes if (mode|>getattr(:r) == origin)))
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
        wannieremptyisometry = globalwannierfunction(extendedemptyseeds,extendedwannierizedemptys[:,resunictcellfockempty])
        wannierfilledisometry = globalwannierfunction(extendedfilledseeds,extendedwannierizedfilleds[:,resunictcellfockfilled])
    
        couriercorrelations = wanniercourierisometry' * correlations * wanniercourierisometry
        couriercorrelationspectrum = couriercorrelations |> crystalspectrum
        purifiedcorrelationspectrum = couriercorrelationspectrum |> roundingpurification
        purifiedcouriercorrelations = purifiedcorrelationspectrum |> CrystalFockMap
    
        
        # emptycorrelations = wannieremptyisometry' * correlations * wannieremptyisometry
        # filledcorrelations = wannierfilledisometry' * correlations * wannierfilledisometry
    
        return Dict(
            :emptyisometry => wannieremptyisometry,
            :filledisometry => wannierfilledisometry,
            :courierisometry => wanniercourierisometry,
            :couriercorrelations => couriercorrelations,
            :correlations => purifiedcouriercorrelations)
    end

    function gmera2ndstep(correlations::FockMap,localcenter1,localcenter2)
        crystal::Crystal = getcrystal(correlations|>getinspace)

        function fulllocalwannierization(correlations::FockMap,localcenter::Offset,selectionstragedy::Function=modeselectionbycount(3))
            shift = [1/2,1/2] ∈ blockedspace
            shiftedunitcell::Region =  Subset(pt+localcenter*2 for pt in (crystal|>getunitcell))
            doubleunitcell::Region = shiftedunitcell+(crystal|>getunitcell)
            shiftedorigunitcell::Region = Subset(pt+localcenter for pt in (blockedcrystal|>getunitcell))
            localregion = intersect(doubleunitcell,shiftedorigunitcell)
            localseedingfock::RegionFock = quantize(localregion,1)
            modebydist = groupmodesbydist(region = localregion,regionfock = localseedingfock,center = localcenter+shift) 
            (emptyseedslist,filledseedslist) = constructfilledandemptyseed(modebydist[length(modebydist)],c6|>recenter(localcenter+shift))
            localiso = localisometries(correlations, localseedingfock, selectionstrategy=selectionstragedy)
    
            emptyseeds = idmap(localseedingfock, localseedingfock)[:,FockSpace(mode for mode in emptyseedslist)] 
            filledseeds = idmap(localseedingfock, localseedingfock)[:,FockSpace(mode for mode in filledseedslist)] 
            courierseeds = idmap(localseedingfock, localseedingfock)[:,sum([FockSpace(mode[2] for mode in modebydist[i]) for i in range(1,length(modebydist)-1)])]

            wannieremptys = _localwannierization(localiso[:empty], emptyseeds)
            wannierfilleds = _localwannierization(localiso[:filled], filledseeds)
            wanniercouriers = _localwannierization(localiso[:courier], courierseeds)
    
            return wannieremptys,wannierfilleds, wanniercouriers,emptyseeds,filledseeds,courierseeds
        end
        
        wannierizedemptys1,wannierizedfilleds1, wannierizedcouriers1,emptyseeds1,filledseeds1,courierseeds1 = fulllocalwannierization(correlations, localcenter1)
        wannierizedemptys2,wannierizedfilleds2, wannierizedcouriers2,emptyseeds2,filledseeds2,courierseeds2 = fulllocalwannierization(correlations, localcenter2)
    
        extendedwannierizedemptys =  wannierizedemptys1 + wannierizedemptys2 
        extendedwannierizedfilleds =  wannierizedfilleds1 + wannierizedfilleds2
        extendedwannierizedcouriers =  wannierizedcouriers1 + wannierizedcouriers2
    
        extendedcourierseeds = courierseeds1 + courierseeds2 
        extendedemptyseeds = emptyseeds1 + emptyseeds2
        extendedfilledseeds = filledseeds1 + filledseeds2 
        
        origin = [0, 0] ∈ blockedspace 
        resunictcellfockcourier = FockSpace(Subset(mode for mode in extendedcourierseeds |> getinspace |> orderedmodes if (mode|>getattr(:r) == origin)))
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
        wannieremptyisometry = globalwannierfunction(extendedemptyseeds,extendedwannierizedemptys[:,resunictcellfockempty])
        wannierfilledisometry = globalwannierfunction(extendedfilledseeds,extendedwannierizedfilleds[:,resunictcellfockfilled])
    
        couriercorrelations = wanniercourierisometry' * correlations * wanniercourierisometry
        couriercorrelationspectrum = couriercorrelations |> crystalspectrum
        purifiedcorrelationspectrum = couriercorrelationspectrum |> roundingpurification
        purifiedcouriercorrelations = purifiedcorrelationspectrum |> CrystalFockMap
    
        
        # emptycorrelations = wannieremptyisometry' * correlations * wannieremptyisometry
        # filledcorrelations = wannierfilledisometry' * correlations * wannierfilledisometry
    
        return Dict(
            :emptyisometry => wannieremptyisometry,
            :filledisometry => wannierfilledisometry,
            :courierisometry => wanniercourierisometry,
            :couriercorrelations => couriercorrelations,
            :correlations => purifiedcouriercorrelations)
    end

    firstcenter1,firstcenter2 = [1/2,0] ∈ blockedspace,[-1/2,0] ∈ blockedspace
    secondcenter1,secondcenter2 = [0,1/2] ∈ blockedspace,[0,-1/2] ∈ blockedspace
    thirdcenter1,thirdcenter2 = [1/2,1/2] ∈ blockedspace,[-1/2,-1/2] ∈ blockedspace
    
    gmera1st = gmera1ststep(blockedcorrelations,firstcenter1,firstcenter2)
    gmera2nd = gmera2ndstep(gmera1st[:correlations],secondcenter1,secondcenter2)
    gmera3rd = gmera2ndstep(gmera2nd[:correlations],thirdcenter1,thirdcenter2)

    return Dict(
            :gmera1stemptyisometry => gmera1st[:emptyisometry],
            :gmera1stfilledisometry => gmera1st[:filledisometry],
            :gmera1stcourierisometry => gmera1st[:courierisometry],
            :gmera1stcorrelations => gmera1st[:correlations],
            :gmera2ndemptyisometry => gmera2nd[:emptyisometry],
            :gmera2ndfilledisometry => gmera2nd[:filledisometry],
            :gmera2ndcourierisometry => gmera2nd[:courierisometry],
            :gmera2ndcorrelations => gmera2nd[:correlations],
            :gmera3rdemptyisometry => gmera3rd[:emptyisometry],
            :gmera3rdfilledisometry => gmera3rd[:filledisometry],
            :gmera3rdcourierisometry => gmera3rd[:courierisometry],
            :gmera3rdcorrelations => gmera3rd[:correlations])
end

rg1 = gmera(correlations,1)

crystalfock = correlations|>getoutspace
scale = Scale([4 0; 0 4], crystalfock|>getcrystal|>getspace)
blocker = @time scale * crystalfock
blockedcorrelations = @time blocker * correlations * blocker'
blockedcrystalfock = blockedcorrelations|>getoutspace
blockedcrystal::Crystal = blockedcrystalfock|>getcrystal
blockedspace::RealSpace = blockedcrystal|>getspace

gmeracorrelations = rg1[:gmera2ndcorrelations]
gmeracrystalfock = gmeracorrelations|>getoutspace
gmeracrystal::Crystal = gmeracrystalfock|>getcrystal
gmeraspace::RealSpace = gmeracrystal|>getspace
thirdcenter1,thirdcenter2 = [1/2,1/2] ∈ blockedspace,[-1/2,-1/2] ∈ blockedspace
shift = [1/2,1/2] ∈ blockedspace

right = [1,0] ∈ blockedspace
up = [0,1] ∈ blockedspace
shifteddigunitcell::Region =  Subset(pt+thirdcenter1*2 for pt in (gmeracrystal|>getunitcell))
shiftedrightunitcell::Region =  Subset(pt+right for pt in (gmeracrystal|>getunitcell))
shiftedupunitcell::Region =  Subset(pt+up for pt in (gmeracrystal|>getunitcell))
qudrpleunitcell::Region = shifteddigunitcell+shiftedrightunitcell+shiftedupunitcell+(gmeracrystal|>getunitcell)
shiftedorigunitcell::Region = Subset(pt+thirdcenter1 for pt in (blockedcrystal|>getunitcell))
localregion = intersect(qudrpleunitcell,shiftedorigunitcell)

localseedingfock::RegionFock = quantize(localregion,1)
modebydist = groupmodesbydist(region = localregion,regionfock = localseedingfock,center = thirdcenter1+shift)
(emptyseedslist,filledseedslist) = constructfilledandemptyseed(modebydist[length(modebydist)],c6|>recenter(thirdcenter1+shift))
localiso = localisometries(gmeracorrelations, localseedingfock, selectionstrategy=modeselectionbycount(3))
emptyseeds = idmap(localseedingfock, localseedingfock)[:,FockSpace(mode for mode in emptyseedslist)] 
filledseeds = idmap(localseedingfock, localseedingfock)[:,FockSpace(mode for mode in filledseedslist)] 
courierseeds = idmap(localseedingfock, localseedingfock)[:,sum([FockSpace(mode[2] for mode in modebydist[i]) for i in range(1,length(modebydist)-1)])]

wannieremptys = _localwannierization(localiso[:empty], emptyseeds)
wannierfilleds = _localwannierization(localiso[:filled], filledseeds)
wanniercouriers = _localwannierization(localiso[:courier], courierseeds)

shiftblockedunitcell = Subset(pt+shift  for pt in (blockedcrystal|>getunitcell))
refdoubleunitcell = shiftblockedunitcell+(blockedcrystal|>getunitcell)
refdoubleunitcell
interest = intersect(refdoubleunitcell,(blockedcrystal|>getunitcell))
shitfedblockedunitcellfock = quantize(shiftblockedunitcell,1)
visualize(regioncorrelations(blockedcorrelations,shitfedblockedunitcellfock)|>eigspech)





blockedunitcellfock = RegionFock(blockedcrystalfock|>unitcellfock)
# visualize(regioncorrelations(blockedcorrelations,blockedunitcellfock)|>eigspech)
localcenter = [1/2,1/2] ∈ blockedspace
localiso = localisometries(blockedcorrelations, blockedunitcellfock, selectionstrategy=modeselectionbycount(3))
modebydist = groupmodesbydist(region = blockedcrystal|>getunitcell,regionfock = blockedunitcellfock,center = localcenter) 
(emptyseedslist,filledseedslist) = constructfilledandemptyseed(modebydist[6],c6|>recenter(localcenter))
emptyseeds = idmap(blockedunitcellfock, blockedunitcellfock)[:,FockSpace(mode for mode in emptyseedslist)] 
filledseeds = idmap(blockedunitcellfock, blockedunitcellfock)[:,FockSpace(mode for mode in filledseedslist)] 
courierseeds = idmap(blockedunitcellfock, blockedunitcellfock)[:,sum([FockSpace(mode[2] for mode in modebydist[i]) for i in range(1,5)])]
wannieremptys = _localwannierization(localiso[:empty], emptyseeds)
wannierfilleds = _localwannierization(localiso[:filled], filledseeds)
wanniercouriers = _localwannierization(localiso[:courier], courierseeds)


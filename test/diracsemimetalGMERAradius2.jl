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


function gmera(correlations::FockMap, step::Int)
    crystalfock = correlations|>getoutspace
    if step == 1
        scale = Scale([6 0; 0 6], crystalfock|>getcrystal|>getspace)
        
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

    function gmera1stor2ndstep(correlations::FockMap,localcenter1,localcenter2,localcenter3,localcenter4,localcenter5)
        crystal::Crystal = getcrystal(correlations|>getinspace)
        space::RealSpace = crystal|>getspace
    
        function fulllocalwannierization(correlations::FockMap,localcenter::Offset,radius::Real=1/3,selectionstragedy::Function=modeselectionbycount(3))
            localregion::Region = getsphericalregionwifcenter(center=localcenter,crystal=crystal, radius=radius, metricspace=space|>orthospace)
            localseedingfock::RegionFock = quantize(localregion,1)
            modebydist = groupmodesbydist(region = localregion,regionfock = localseedingfock,center = localcenter) 
            (emptyseedslist,filledseedslist) = constructfilledandemptyseed(modebydist[3],c6|>recenter(localcenter))
            localiso = localisometries(correlations, localseedingfock, selectionstrategy=selectionstragedy)
    
            emptyseeds = idmap(localseedingfock, localseedingfock)[:,FockSpace(mode for mode in emptyseedslist)] 
            filledseeds = idmap(localseedingfock, localseedingfock)[:,FockSpace(mode for mode in filledseedslist)] 
            courierseeds = idmap(localseedingfock, localseedingfock)[:,FockSpace(mode[2] for mode in modebydist[1])+FockSpace(mode[2] for mode in modebydist[2])]

            wannieremptys = _localwannierization(localiso[:empty], emptyseeds)
            wannierfilleds = _localwannierization(localiso[:filled], filledseeds)
            wanniercouriers = _localwannierization(localiso[:courier], courierseeds)
    
            return wannieremptys,wannierfilleds, wanniercouriers,emptyseeds,filledseeds,courierseeds
        end
        
        wannierizedemptys1,wannierizedfilleds1, wannierizedcouriers1,emptyseeds1,filledseeds1,courierseeds1 = fulllocalwannierization(correlations, localcenter1)
        wannierizedemptys2,wannierizedfilleds2, wannierizedcouriers2,emptyseeds2,filledseeds2,courierseeds2 = fulllocalwannierization(correlations, localcenter2)
        wannierizedemptys3,wannierizedfilleds3, wannierizedcouriers3,emptyseeds3,filledseeds3,courierseeds3 = fulllocalwannierization(correlations, localcenter3)
        wannierizedemptys4,wannierizedfilleds4, wannierizedcouriers4,emptyseeds4,filledseeds4,courierseeds4 = fulllocalwannierization(correlations, localcenter4)
        wannierizedemptys5,wannierizedfilleds5, wannierizedcouriers5,emptyseeds5,filledseeds5,courierseeds5 = fulllocalwannierization(correlations, localcenter5)
    
        extendedwannierizedemptys =  wannierizedemptys1 + wannierizedemptys2 + wannierizedemptys3 + wannierizedemptys4 + wannierizedemptys5
        extendedwannierizedfilleds =  wannierizedfilleds1 + wannierizedfilleds2 + wannierizedfilleds3 + wannierizedfilleds4 + wannierizedfilleds5
        extendedwannierizedcouriers =  wannierizedcouriers1 + wannierizedcouriers2 + wannierizedcouriers3 + wannierizedcouriers4 + wannierizedcouriers5
    
        extendedcourierseeds = courierseeds1 + courierseeds2 + courierseeds3 + courierseeds4 + courierseeds5
        extendedemptyseeds = emptyseeds1 + emptyseeds2 + emptyseeds3 + emptyseeds4 + emptyseeds5
        extendedfilledseeds = filledseeds1 + filledseeds2 + filledseeds3 + filledseeds4 + filledseeds5
        
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

    function gmera3rdstep(correlations::FockMap,localcenter1,localcenter2,localcenter3,localcenter4,localcenter5,localcenter6)
        crystal::Crystal = getcrystal(correlations|>getinspace)
        space::RealSpace = crystal|>getspace
    
        function fulllocalwannierization(correlations::FockMap,localcenter::Offset,radius::Real=1/3,selectionstragedy::Function=modeselectionbycount(3))
            localregion::Region = getsphericalregionwifcenter(center=localcenter,crystal=crystal, radius=radius, metricspace=space|>orthospace)
            localseedingfock::RegionFock = quantize(localregion,1)
            modebydist = groupmodesbydist(region = localregion,regionfock = localseedingfock,center = localcenter) 
            (emptyseedslist,filledseedslist) = constructfilledandemptyseed(modebydist[2],c6|>recenter(localcenter))
            localiso = localisometries(correlations, localseedingfock, selectionstrategy=selectionstragedy)
    
            emptyseeds = idmap(localseedingfock, localseedingfock)[:,FockSpace(mode for mode in emptyseedslist)] 
            filledseeds = idmap(localseedingfock, localseedingfock)[:,FockSpace(mode for mode in filledseedslist)] 
    
            courierseeds = idmap(localseedingfock, localseedingfock)[:,FockSpace(mode[2] for mode in modebydist[1])]
            wannieremptys = _localwannierization(localiso[:empty], emptyseeds)
            wannierfilleds = _localwannierization(localiso[:filled], filledseeds)
            wanniercouriers = _localwannierization(localiso[:courier], courierseeds)
    
            return wannieremptys,wannierfilleds, wanniercouriers,emptyseeds,filledseeds,courierseeds
        end
        
        wannierizedemptys1,wannierizedfilleds1, wannierizedcouriers1,emptyseeds1,filledseeds1,courierseeds1 = fulllocalwannierization(correlations, localcenter1)
        wannierizedemptys2,wannierizedfilleds2, wannierizedcouriers2,emptyseeds2,filledseeds2,courierseeds2 = fulllocalwannierization(correlations, localcenter2)
        wannierizedemptys3,wannierizedfilleds3, wannierizedcouriers3,emptyseeds3,filledseeds3,courierseeds3 = fulllocalwannierization(correlations, localcenter3)
        wannierizedemptys4,wannierizedfilleds4, wannierizedcouriers4,emptyseeds4,filledseeds4,courierseeds4 = fulllocalwannierization(correlations, localcenter4)
        wannierizedemptys5,wannierizedfilleds5, wannierizedcouriers5,emptyseeds5,filledseeds5,courierseeds5 = fulllocalwannierization(correlations, localcenter5)
        wannierizedemptys6,wannierizedfilleds6, wannierizedcouriers6,emptyseeds6,filledseeds6,courierseeds6 = fulllocalwannierization(correlations, localcenter6)
    
        extendedwannierizedemptys =  wannierizedemptys1 + wannierizedemptys2 + wannierizedemptys3 + wannierizedemptys4 + wannierizedemptys5 + wannierizedemptys6
        extendedwannierizedfilleds =  wannierizedfilleds1 + wannierizedfilleds2 + wannierizedfilleds3 + wannierizedfilleds4 + wannierizedfilleds5 + wannierizedfilleds6
        extendedwannierizedcouriers =  wannierizedcouriers1 + wannierizedcouriers2 + wannierizedcouriers3 + wannierizedcouriers4 + wannierizedcouriers5 + wannierizedcouriers6
    
        extendedcourierseeds = courierseeds1 + courierseeds2 + courierseeds3 + courierseeds4 + courierseeds5 + courierseeds6
        extendedemptyseeds = emptyseeds1 + emptyseeds2 + emptyseeds3 + emptyseeds4 + emptyseeds5 + emptyseeds6
        extendedfilledseeds = filledseeds1 + filledseeds2 + filledseeds3 + filledseeds4 + filledseeds5 + filledseeds6
        
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

    firstcenter1,firstcenter2,firstcenter3,firstcenter4,firstcenter5 = [1/3,0] ∈ blockedspace,[0,1/3] ∈ blockedspace,[1/3,1] ∈ blockedspace,[2/3,2/3] ∈ blockedspace,[1,1/3] ∈ blockedspace
    secondcenter1,secondcenter2,secondcenter3,secondcenter4,secondcenter5 = [2/3,0] ∈ blockedspace,[0,2/3] ∈ blockedspace,[1/3,1/3] ∈ blockedspace,[1,2/3] ∈ blockedspace,[2/3,1] ∈ blockedspace
    thirdcenter1,thirdcenter2,thirdcenter3,thirdcenter4,thirdcenter5,thirdcenter6 = [0,0] ∈ blockedspace,[0,1] ∈ blockedspace,[1,0] ∈ blockedspace,[1/3,2/3] ∈ blockedspace,[2/3,1/3] ∈ blockedspace,[1, 1] ∈ blockedspace
    
    gmera1ststep = gmera1stor2ndstep(blockedcorrelations,firstcenter1,firstcenter2,firstcenter3,firstcenter4,firstcenter5)
    gmera2ndstep = gmera1stor2ndstep(gmera1ststep[:correlations],secondcenter1,secondcenter2,secondcenter3,secondcenter4,secondcenter5)
    gmera3rdstep = gmera3rdstep(gmera2ndstep[:correlations],thirdcenter1,thirdcenter2,thirdcenter3,thirdcenter4,thirdcenter5,thirdcenter6)

    return Dict(
            :gmera1stemptyisometry => gmera1ststep[:emptyisometry],
            :gmera1stfilledisometry => gmera1ststep[:filledisometry],
            :gmera1stcourierisometry => gmera1ststep[:courierisometry],
            :gmera1stcorrelations => gmera1ststep[:correlations],
            :gmera2ndemptyisometry => gmera2ndstep[:emptyisometry],
            :gmera2ndfilledisometry => gmera2ndstep[:filledisometry],
            :gmera2ndcourierisometry => gmera2ndstep[:courierisometry],
            :gmera2ndcorrelations => gmera2ndstep[:correlations],
            :gmera3rdemptyisometry => gmera3rdstep[:emptyisometry],
            :gmera3rdfilledisometry => gmera3rdstep[:filledisometry],
            :gmera3rdcourierisometry => gmera3rdstep[:courierisometry],
            :gmera3rdcorrelations => gmera3rdstep[:correlations])
end

rg1 = gmera(correlations,1)
rg2 = gmera(rg1[:gmera3rdcorrelations],2)
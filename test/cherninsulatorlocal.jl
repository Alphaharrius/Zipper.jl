using Plotly, SmithNormalForm, LinearAlgebra, OrderedCollections, SparseArrays, Combinatorics
using Revise
using Zipper
using DataFrames

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

function findlocalseeds(;
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

    localcorrelations::FockMap = statecorrelations[regionfock, regionfock]
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

function _findlocalseeds(;
    statecorrelations::FockMap, regionmodes::Subset{Mode},
    spectrumextractpredicate::Function = v -> v < 1e-2)
 
    localcorrelations::FockMap = statecorrelations[regionmodes |> FockSpace, regionmodes|> FockSpace] 
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

function circularfilter(mode::Mode, center::Offset, radius::Real = 1.5)::Bool
    return norm(((mode-center) |> getpos |> euclidean))<=radius
end

function circularregionmodes(center::Offset, physicalmodes::Subset{Mode}, radius::Number, crystal:: Crystal)::Subset{Mode}
    currentspace::RealSpace = correlations.outspace |> getcrystal |> getspace |> orthogonalspace
    crystalsize = crystal |> size
    # periodicnorm =  m -> lineartransform(currentspace, takeperiodic(m |> getpos,crystalsize)) |> norm
     
    # function takeperiodic(p, crystalsize)
    #     return (p |> getspace)&[min(vectno,24-vectno) for (size,vectno) in zip(crystalsize, vec(p))]
    # end
    return filter(p -> circularfilter(p, center,radius), physicalmodes)
end


function localmodesgrouping(localcorrelation::FockMap, threshold::Float64)::Dict{Symbol, FockMap} 
    spectrum::EigenSpectrum = localcorrelation |> eigspech
    filledmodes::Subset{Mode} = Subset(emode for (emode, eval) in (filter(p -> p.second < threshold, spectrum |> geteigenvalues)))
    emptymodes::Subset{Mode} = Subset(emode for (emode, eval) in filter(p -> p.second > 1.0 - threshold, spectrum |> geteigenvalues))
    frozenmodes::Subset{Mode} = Subset(emode for (emode, eval) in filter(p -> (p.second > 1.0 - threshold || p.second < threshold), spectrum |> geteigenvalues))
    couriermodes::Subset{Mode} = Subset(emode for (emode, eval) in filter(p -> threshold < p.second < 1.0 - threshold, spectrum |> geteigenvalues))
    return Dict(:filled => columns(spectrum |> geteigenvectors, FockSpace(filledmodes)), :empty => columns(spectrum |> geteigenvectors, FockSpace(emptymodes)),
    :frozen => columns(spectrum |> geteigenvectors, FockSpace(frozenmodes)), :courier => columns(spectrum |> geteigenvectors, FockSpace(couriermodes)))
end

function localregioninspection(center::Offset, physicalmodes::Subset{Mode}, radius::Number, crystal:: Crystal)::Tuple{Subset{Offset},FockSpace}
    seedingmodes::Subset{Mode} = circularregionmodes(center, physicalmodes, radius, crystal)
    seedingregion::Subset{Offset} = Subset(m |> getpos for m in seedingmodes)
    seedingfock::FockSpace = RegionFock(seedingmodes)
    return seedingregion,seedingfock
end

function localwannierization(localbasis::FockMap, localseeds::FockMap, svdorthothreshold::Number = 1e-1)::FockMap
    U, Σ, Vt = svd(localbasis'*localseeds)
    minsvdvalue::Number = minimum(v for (_, v) in Σ)
    println("min svdvalue", minsvdvalue)
    precarioussvdvalues::Vector = []
    if minsvdvalue < svdorthothreshold
        push!(precarioussvdvalues, minsvdvalue)
    end
    unitary::FockMap = U * Vt
    wannierizedbasis = localbasis*unitary
    wannierizedbasis = wannierizedbasis*spatialmap(wannierizedbasis)'
    return wannierizedbasis
end

function _spatialmap(fockmap::FockMap)::FockMap
    function _spatialinmode(colmap::FockMap, ind::Integer)
        inmode::Mode = colmap |> getinspace |> first
        absmap::FockMap = colmap |> abs 
        modecenter::Offset = sort(absmap |> Zipper.columnspec, by=p->p.second |> real) |> last |> first |> getpos
        basis::Offset = modecenter |> basispoint
        offset::Offset = modecenter - basis
        return inmode |> setattr(:r => offset) |> setattr(:b => basis) |> setattr(:ind => ind)
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

function groupmodesbydist(;
    region::Subset{Point{RealSpace}},
    regionfock::FockSpace,
    center::Point,
    samedistancethreshold::Int = 8)

    visualspace = region |> getspace |> euclidean
    distancewithmode = sort(((norm(lineartransform(visualspace, (mode |> getpos)-center) |> vec),mode) for mode in regionfock), by = first, rev = true)
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


∈(vec, space::AffineSpace)::Point = Point(vec |> collect, space)

triangular = RealSpace([sqrt(3)/2 -1/2; 0. 1.]')
kspace = convert(MomentumSpace, triangular)

pa = [1/3, 2/3] ∈ triangular
pb = [2/3, 1/3] ∈ triangular
pc = (pa + pb) / 2
spatialsnappingcalibration((pa, pb, pc))

c6 = pointgrouptransform([cos(π/3) -sin(π/3); sin(π/3) cos(π/3)])

unitcell = Subset(pa, pb)
crystal = Crystal(unitcell, [96, 96])
reciprocalhashcalibration(crystal.sizes)

modes::Subset{Mode} = quantize(:b, unitcell, 1)
m0, m1 = members(modes)

tₙ = ComplexF64(-1)
tₕ = 0.1im

nearestneighbor = [
    (m0, m1) => tₙ,
    (m0, setattr(m1, :r => Point([-1, 0], triangular))) => tₙ,
    (m0, setattr(m1, :r => Point([0, 1], triangular))) => tₙ]

haldane = [
    (m0, setattr(m0, :r => Point([1, 1], triangular))) => tₕ,
    (m0, setattr(m0, :r => Point([-1, 0], triangular))) => tₕ,
    (m0, setattr(m0, :r => Point([0, -1], triangular))) => tₕ,
    (m1, setattr(m1, :r => Point([1, 1], triangular))) => -tₕ,
    (m1, setattr(m1, :r => Point([-1, 0], triangular))) => -tₕ,
    (m1, setattr(m1, :r => Point([0, -1], triangular))) => -tₕ]

bonds::FockMap = bondmap([nearestneighbor..., haldane...])

energyspectrum = computeenergyspectrum(bonds, crystal=crystal)
energyspectrum |> visualize

groundstates::CrystalSpectrum = groundstatespectrum(energyspectrum, perunitcellfillings=1)
groundstates |> visualize
groundstateprojector = groundstates |> crystalprojector

C = idmap(groundstateprojector.outspace) - groundstateprojector

C |> crystalspectrum |> visualize

correlations = C

crystalfock = correlations.outspace

scale = Scale([2 0; 0 2], crystalfock |> getcrystal |> getspace)
blockresult = blocking(:action => scale, :correlations => correlations, :crystal => crystalfock |> getcrystal)

blockedcrystal::Crystal = blockresult[:crystal]
blockedcorrelations::FockMap = blockresult[:correlations]

crystalpoints::Subset{Offset} = latticepoints(blockedcrystal)
blockedmodes::Subset{Mode} = quantize(:b, blockedcrystal.unitcell, 1)
physicalmodes::Subset{Mode} = spanoffset(blockedmodes, crystalpoints)
scaledtriangular = scale*triangular 

center_pt = [-2,0]
center = Point(center_pt, scaledtriangular)

localregion,localfock = localregioninspection(center , physicalmodes, 2,blockedcrystal)
visualize(localregion)

localcorrelation = regioncorrelations(blockedcorrelations,localfock)
localeigspec = localcorrelation |> eigspech 
visualize(localeigspec)

firstlayer = Set([[-2,-2], [-1,0], [0,2], [0,-1], [1,1], [2,0]])
secondlayer = Set([[-1,-1], [0,1], [1,0]])
thirdlayer = Set([[0,0]])

firstlayers = Set([[-2,-2], [-1,0], [0,2], [0,-1], [1,1], [2,0]])
secondlayers = Set([[-1,-1], [0,1], [1,0]])
thirdlayers = Set([[0,0]])

for i in range(-5,5)
    for j in range(-5,5)
        firstlayers = union(firstlayers,Set([r+i*[1,2]+j*[2,1] for r in firstlayer]))
        secondlayers = union(secondlayers,Set([r+i*[1,2]+j*[2,1] for r in secondlayer]))
        thirdlayers = union(thirdlayers,Set([r+i*[1,2]+j*[2,1] for r in thirdlayer]))
    end
end

function locaclRG(center_pt)::Tuple{FockMap,FockMap}
    center = Point(center_pt, scaledtriangular)

    localregion,localfock = localregioninspection(center , physicalmodes, 2,blockedcrystal)
    visualize(localregion)

    localcorrelation = regioncorrelations(blockedcorrelations,localfock)
    localeigspec = localcorrelation |> eigspech 
    visualize(localeigspec)

    localmodesdict = localmodesgrouping(localcorrelation, 0.001)
    modebydist = groupmodesbydist(region = localregion,regionfock = localfock,center = center)  

    # finding seeds for local frozen (6 modes at the center)
    frozenseedingfocks = (frozenseedingfock for (_,frozenseedingfock) in (localregioninspection(modewithdist[2] |> getpos, physicalmodes, 0.1, blockedcrystal) for modewithdist in modebydist[3]))
    frozenseeds = sum([reduce(+ ,FockMap(v,inspace = v |> getinspace |> orderedmodes |> setattr(:dumind => n) |> FockSpace,performpermute = false) for (n, (_ , v)) in enumerate(_findlocalspstates(statecorrelations = blockedcorrelations, regionfock = frozenseedingfock, 
    spectrumextractpredicate = v -> true, symmetry = identitytransform(2)))) for frozenseedingfock in frozenseedingfocks])

    extendedfrozenseeds = idmap(localfock, localfock)[:,frozenseeds |> getoutspace] * frozenseeds

    wannierizedfrozens = _localwannierization(localmodesdict[:frozen], extendedfrozenseeds)

    # finding seeds for local courier (18 (12+6) modes at "boundary")
    courierseedingfocks1 = (courierseedingfock for (_,courierseedingfock) in (localregioninspection(modewithdist[2] |> getpos, physicalmodes, 0.1, blockedcrystal) for modewithdist in modebydist[1]))
    courierseedingfocks2 = (courierseedingfock for (_,courierseedingfock) in (localregioninspection(modewithdist[2] |> getpos, physicalmodes, 0.1, blockedcrystal) for modewithdist in modebydist[2]))

    courierseeds1 = sum([reduce(+ ,FockMap(v,inspace = v |> getinspace |> orderedmodes |> setattr(:dumind => n) |> FockSpace,performpermute = false) for (n, (_ , v)) in enumerate(_findlocalspstates(statecorrelations = blockedcorrelations, regionfock = courierseedingfock, 
    spectrumextractpredicate = v -> true, symmetry = identitytransform(2)))) for courierseedingfock in courierseedingfocks1])
    courierseeds2 = sum([reduce(+ ,FockMap(v,inspace = v |> getinspace |> orderedmodes |> setattr(:dumind => n) |> FockSpace,performpermute = false) for (n, (_ , v)) in enumerate(_findlocalspstates(statecorrelations = blockedcorrelations, regionfock = courierseedingfock, 
    spectrumextractpredicate = v -> true, symmetry = identitytransform(2)))) for courierseedingfock in courierseedingfocks2])

    courierseeds = courierseeds1 + courierseeds2
    extendedcourierseeds = idmap(localfock, localfock)[:, courierseeds |> getoutspace] * courierseeds

    wannierizedcouriers = _localwannierization(localmodesdict[:courier], extendedcourierseeds)

    return wannierizedfrozens, wannierizedcouriers
end

function locaclRGsecond(center_pt)::Tuple{FockMap,FockMap}
    center = Point(center_pt, scaledtriangular)
    # inspecting local region spectrum
    transformedphysicalmodes = couriercorrelation |> getoutspace |> orderedmodes
    
    trsasnformedlocalregion, trsasnformedlocalfock = localregioninspection(center ,transformedphysicalmodes, 2, blockedcrystal)
    visualize(trsasnformedlocalregion, title="l", visualspace=euclidean(RealSpace, 2))
    
    transformedlocaleigspec = couriercorrelation[trsasnformedlocalfock, trsasnformedlocalfock] |> eigspech 
    visualize(transformedlocaleigspec)
    
    
    transformedlocalmodesdict = localmodesgrouping(couriercorrelation[trsasnformedlocalfock, trsasnformedlocalfock], 0.01)
    transformedmodebydist = groupmodesbydist(region = trsasnformedlocalregion,regionfock = trsasnformedlocalfock,center = center)
    
    # finding seeds for local frozen (6 modes at the center)
    frozenseedingfocksRG = (frozenseedingfock for (_,frozenseedingfock) in (localregioninspection(modewithdist[2] |> getpos, transformedphysicalmodes, 0.1, blockedcrystal) for modewithdist in transformedmodebydist[3]))
    frozenseedsRG = sum([reduce(+ ,FockMap(v,inspace = v |> getinspace |> orderedmodes |> setattr(:dumind => n) |> FockSpace,performpermute = false) for (n, (_ , v)) in enumerate(findlocalseeds(statecorrelations = couriercorrelation, regionfock = frozenseedingfockRG, 
    spectrumextractpredicate = v -> true, symmetry = identitytransform(2)))) for frozenseedingfockRG in frozenseedingfocksRG])

    extendedfrozenseedsRG = idmap(trsasnformedlocalfock, trsasnformedlocalfock)[:,frozenseedsRG |> getoutspace] * frozenseedsRG

    wannierizedfrozensRG = _localwannierization(transformedlocalmodesdict[:frozen], extendedfrozenseedsRG)

    # finding seeds for local courier (18 (12+6) modes at "boundary")
    courierseedingfocks1RG = (courierseedingfock for (_,courierseedingfock) in (localregioninspection(modewithdist[2] |> getpos, transformedphysicalmodes, 0.1, blockedcrystal) for modewithdist in transformedmodebydist[1]))
    courierseedingfocks2RG = (courierseedingfock for (_,courierseedingfock) in (localregioninspection(modewithdist[2] |> getpos, transformedphysicalmodes, 0.1, blockedcrystal) for modewithdist in transformedmodebydist[2]))

    courierseeds1RG = sum([reduce(+ ,FockMap(v,inspace = v |> getinspace |> orderedmodes |> setattr(:dumind => n) |> FockSpace,performpermute = false) for (n, (_ , v)) in enumerate(findlocalseeds(statecorrelations = couriercorrelation, regionfock = courierseedingfockRG, 
    spectrumextractpredicate = v -> true, symmetry = identitytransform(2)))) for courierseedingfockRG in courierseedingfocks1RG])
    courierseeds2RG = sum([reduce(+ ,FockMap(v,inspace = v |> getinspace |> orderedmodes |> setattr(:dumind => n) |> FockSpace,performpermute = false) for (n, (_ , v)) in enumerate(findlocalseeds(statecorrelations = couriercorrelation, regionfock = courierseedingfockRG, 
    spectrumextractpredicate = v -> true, symmetry = identitytransform(2)))) for courierseedingfockRG in courierseedingfocks2RG])

    courierseedsRG = courierseeds1RG + courierseeds2RG
    extendedcourierseedsRG = idmap(trsasnformedlocalfock, trsasnformedlocalfock)[:, courierseedsRG |> getoutspace] * courierseedsRG
    

    # Wannierization for couriers
    wannierizedcouriersRG = _localwannierization(transformedlocalmodesdict[:courier], extendedcourierseedsRG)
    
    return wannierizedfrozensRG, wannierizedcouriersRG
end

function locaclRGthird(center_pt)::Tuple{FockMap,FockMap}
    center = Point(center_pt, scaledtriangular)
    # inspecting local region spectrum
    transformedphysicalmodessecond = couriercorrelationRG |> getoutspace |> orderedmodes
    
    trsasnformedlocalregionsecond, trsasnformedlocalfocksecond = localregioninspection(center,transformedphysicalmodessecond, 2, blockedcrystal)
    visualize(trsasnformedlocalregionsecond, title="l", visualspace=euclidean(RealSpace, 2))
    
    transformedlocaleigspec = couriercorrelationRG[trsasnformedlocalfocksecond, trsasnformedlocalfocksecond] |> eigspech 
    visualize(transformedlocaleigspec)
    
    transformedlocalmodesdictsecond = localmodesgrouping(couriercorrelationRG[trsasnformedlocalfocksecond, trsasnformedlocalfocksecond], 0.002)
    transformedmodebydistsecond = groupmodesbydist(region = trsasnformedlocalregionsecond,regionfock = trsasnformedlocalfocksecond,center = center)
    
    # finding seeds for local frozen (6 modes at the center)
    frozenseedingfocksRGsecond = (frozenseedingfock for (_,frozenseedingfock) in (localregioninspection(modewithdist[2] |> getpos, transformedphysicalmodessecond, 0.1, blockedcrystal) for modewithdist in transformedmodebydistsecond[2]))
    frozenseedsRGsecond = sum([reduce(+ ,FockMap(v,inspace = v |> getinspace |> orderedmodes |> setattr(:dumind => n) |> FockSpace,performpermute = false) for (n, (_ , v)) in enumerate(findlocalseeds(statecorrelations = couriercorrelationRG, regionfock = frozenseedingfockRG, 
    spectrumextractpredicate = v -> true, symmetry = identitytransform(2)))) for frozenseedingfockRG in frozenseedingfocksRGsecond])

    extendedfrozenseedsRGsecond = idmap(trsasnformedlocalfocksecond, trsasnformedlocalfocksecond)[:,frozenseedsRGsecond |> getoutspace] * frozenseedsRGsecond

    wannierizedfrozensRG = _localwannierization(transformedlocalmodesdictsecond[:frozen], extendedfrozenseedsRGsecond)

    # finding seeds for local courier (18 (12+6) modes at "boundary")
    courierseedingfocksRGsecond = (courierseedingfock for (_,courierseedingfock) in (localregioninspection(modewithdist[2] |> getpos, transformedphysicalmodessecond, 0.1, blockedcrystal) for modewithdist in transformedmodebydistsecond[1]))
    
    courierseedsRGsecond = sum([reduce(+ ,FockMap(v,inspace = v |> getinspace |> orderedmodes |> setattr(:dumind => n) |> FockSpace,performpermute = false) for (n, (_ , v)) in enumerate(findlocalseeds(statecorrelations = couriercorrelationRG, regionfock = courierseedingfockRG, 
    spectrumextractpredicate = v -> true, symmetry = identitytransform(2)))) for courierseedingfockRG in courierseedingfocksRGsecond])

    extendedcourierseedsRGsecond = idmap(trsasnformedlocalfocksecond, trsasnformedlocalfocksecond)[:, courierseedsRGsecond |> getoutspace] * courierseedsRGsecond
    

    # Wannierization for couriers
    wannierizedcouriersRG = _localwannierization(transformedlocalmodesdictsecond[:courier], extendedcourierseedsRGsecond)
    
    return wannierizedfrozensRG, wannierizedcouriersRG
end

wannier_fns = [locaclRG(center_pt) for center_pt in firstlayers]
localunitarylist = [wannierizedfrozen+wannierizedcourier for (wannierizedfrozen, wannierizedcourier) in wannier_fns]
localunitary = sum(localunitarylist)
wanniercourierlist = [wannierizedcourier for (_, wannierizedcourier) in wannier_fns]
wanniercouriers = sum(wanniercourierlist)

couriercorrelation = wanniercouriers'*regioncorrelations(blockedcorrelations,localunitary |> getoutspace)*wanniercouriers

wannier_fnsRG = [locaclRGsecond(center_pt) for center_pt in secondlayers]
localunitarylistRG = [wannierizedfrozen+wannierizedcourier for (wannierizedfrozen, wannierizedcourier) in wannier_fnsRG]
localunitaryRG = sum(localunitarylistRG)
wanniercourierlistRG = [wannierizedcourier for (_, wannierizedcourier) in wannier_fnsRG]
wanniercouriersRG = sum(wanniercourierlistRG)

couriercorrelationRG = wanniercouriersRG'*couriercorrelation[localunitaryRG |> getoutspace, localunitaryRG |> getoutspace]*wanniercouriersRG

wannier_fnsRGRG = [locaclRGthird(center_pt) for center_pt in thirdlayers]
localunitarylistRGRG = [wannierizedfrozen+wannierizedcourier for (wannierizedfrozen, wannierizedcourier) in wannier_fnsRGRG]
localunitaryRGRG = sum(localunitarylistRGRG)
wanniercourierlistRGRG = [wannierizedcourier for (_, wannierizedcourier) in wannier_fnsRGRG]
wanniercouriersRGRG = sum(wanniercourierlistRGRG)

couriercorrelationRGRG = wanniercouriersRGRG'*couriercorrelationRG[localunitaryRGRG |> getoutspace, localunitaryRGRG |> getoutspace]*wanniercouriersRGRG

# inspecting local region spectrum
center_pt = [0,0]
center = Point(center_pt, scaledtriangular)
transformedphysicalmodesnext = couriercorrelationRGRG |> getoutspace |> orderedmodes

trsasnformedlocalregionnext, trsasnformedlocalfocknext = localregioninspection(center,transformedphysicalmodesnext, 4, blockedcrystal)
visualize(trsasnformedlocalregionnext, title="l", visualspace=euclidean(RealSpace, 2))

transformedlocaleigspec = couriercorrelationRGRG[trsasnformedlocalfocknext, trsasnformedlocalfocknext] |> eigspech 
visualize(transformedlocaleigspec)


firstlayernext = Set([[-4,-4], [-2,0], [0,4], [0,-2], [2,2], [4,0]])
secondlayernext = Set([[-2,-2], [0,2], [2,0]])
thirdlayernext = Set([[0,0]])

firstlayersnext = Set([[-4,-4], [-2,0], [0,4], [0,-2], [2,2], [4,0]])
secondlayersnext = Set([[-2,-2], [0,2], [2,0]])
thirdlayersnext = Set([[0,0]])

# for i in range(-1,1)
#     for j in range(-1,1)
#         firstlayersnext = union(firstlayersnext,Set([r+i*[1,2]+j*[2,1] for r in firstlayer]))
#         secondlayersnext = union(secondlayersnext,Set([r+i*[1,2]+j*[2,1] for r in secondlayer]))
#         thirdlayersnext = union(thirdlayersnext,Set([r+i*[1,2]+j*[2,1] for r in thirdlayer]))
#     end
# end

function locaclRGnext(center_pt)::Tuple{FockMap,FockMap}
    center = Point(center_pt, scaledtriangular)

    physicalmodesnext = couriercorrelationRGRG |> getoutspace |> orderedmodes

    localregionnext, localfocknext = localregioninspection(center,physicalmodesnext, 4, blockedcrystal)
    visualize(localregionnext, title="l", visualspace=euclidean(RealSpace, 2))

    localeigspec = couriercorrelationRGRG[localfocknext, localfocknext] |> eigspech 
    visualize(localeigspec)

    localmodesdict = localmodesgrouping(couriercorrelationRGRG[localfocknext, localfocknext], 0.002)
    modebydistnext = groupmodesbydist(region = localregionnext,regionfock = localfocknext,center = center)  

    # finding seeds for local frozen (6 modes at the center)
    frozenseedingfocks = (frozenseedingfock for (_,frozenseedingfock) in (localregioninspection(modewithdist[2] |> getpos, physicalmodesnext, 0.1, blockedcrystal) for modewithdist in modebydistnext[3]))
    frozenseeds = sum([reduce(+ ,FockMap(v,inspace = v |> getinspace |> orderedmodes |> setattr(:dumind => n) |> FockSpace,performpermute = false) for (n, (_ , v)) in enumerate(findlocalseeds(statecorrelations = couriercorrelationRGRG, regionfock = frozenseedingfock, 
    spectrumextractpredicate = v -> true, symmetry = identitytransform(2)))) for frozenseedingfock in frozenseedingfocks])

    extendedfrozenseeds = idmap(localfocknext, localfocknext)[:,frozenseeds |> getoutspace] * frozenseeds

    wannierizedfrozens = _localwannierization(localmodesdict[:frozen], extendedfrozenseeds)

    # finding seeds for local courier (18 (12+6) modes at "boundary")
    courierseedingfocks1 = (courierseedingfock for (_,courierseedingfock) in (localregioninspection(modewithdist[2] |> getpos, physicalmodesnext, 0.1, blockedcrystal) for modewithdist in modebydistnext[1]))
    courierseedingfocks2 = (courierseedingfock for (_,courierseedingfock) in (localregioninspection(modewithdist[2] |> getpos, physicalmodesnext, 0.1, blockedcrystal) for modewithdist in modebydistnext[2]))

    courierseeds1 = sum([reduce(+ ,FockMap(v,inspace = v |> getinspace |> orderedmodes |> setattr(:dumind => n) |> FockSpace,performpermute = false) for (n, (_ , v)) in enumerate(findlocalseeds(statecorrelations = couriercorrelationRGRG, regionfock = courierseedingfock, 
    spectrumextractpredicate = v -> true, symmetry = identitytransform(2)))) for courierseedingfock in courierseedingfocks1])
    courierseeds2 = sum([reduce(+ ,FockMap(v,inspace = v |> getinspace |> orderedmodes |> setattr(:dumind => n) |> FockSpace,performpermute = false) for (n, (_ , v)) in enumerate(findlocalseeds(statecorrelations = couriercorrelationRGRG, regionfock = courierseedingfock, 
    spectrumextractpredicate = v -> true, symmetry = identitytransform(2)))) for courierseedingfock in courierseedingfocks2])

    courierseeds = courierseeds1 + courierseeds2
    extendedcourierseeds = idmap(localfocknext, localfocknext)[:, courierseeds |> getoutspace] * courierseeds

    wannierizedcouriers = _localwannierization(localmodesdict[:courier], extendedcourierseeds)

    return wannierizedfrozens, wannierizedcouriers
end

function locaclRGsecondnext(center_pt)::Tuple{FockMap,FockMap}
    center = Point(center_pt, scaledtriangular)
    # inspecting local region spectrum
    transformedphysicalmodesnext = couriercorrelationnext |> getoutspace |> orderedmodes
    
    transformedlocalregionnext, transformedlocalfocknext = localregioninspection(center ,transformedphysicalmodesnext, 4, blockedcrystal)
    visualize(transformedlocalregionnext, title="l", visualspace=euclidean(RealSpace, 2))
    
    transformedlocaleigspecnext = couriercorrelationnext[transformedlocalfocknext, transformedlocalfocknext] |> eigspech 
    visualize(transformedlocaleigspecnext)

    transformedlocalmodesdict = localmodesgrouping(couriercorrelationnext[transformedlocalfocknext, transformedlocalfocknext], 0.003)
    transformedmodebydist = groupmodesbydist(region =  transformedlocalregionnext,regionfock = transformedlocalfocknext,center = center)
    
    # finding seeds for local frozen (6 modes at the center)
    frozenseedingfocksRG = (frozenseedingfock for (_,frozenseedingfock) in (localregioninspection(modewithdist[2] |> getpos, transformedphysicalmodesnext, 0.1, blockedcrystal) for modewithdist in transformedmodebydist[3]))
    frozenseedsRG = sum([reduce(+ ,FockMap(v,inspace = v |> getinspace |> orderedmodes |> setattr(:dumind => n) |> FockSpace,performpermute = false) for (n, (_ , v)) in enumerate(findlocalseeds(statecorrelations = couriercorrelationnext, regionfock = frozenseedingfockRG, 
    spectrumextractpredicate = v -> true, symmetry = identitytransform(2)))) for frozenseedingfockRG in frozenseedingfocksRG])

    extendedfrozenseedsRG = idmap(transformedlocalfocknext, transformedlocalfocknext)[:,frozenseedsRG |> getoutspace] * frozenseedsRG

    wannierizedfrozensRG = _localwannierization(transformedlocalmodesdict[:frozen], extendedfrozenseedsRG)

    # finding seeds for local courier (18 (12+6) modes at "boundary")
    courierseedingfocks1RG = (courierseedingfock for (_,courierseedingfock) in (localregioninspection(modewithdist[2] |> getpos, transformedphysicalmodesnext, 0.1, blockedcrystal) for modewithdist in transformedmodebydist[1]))
    courierseedingfocks2RG = (courierseedingfock for (_,courierseedingfock) in (localregioninspection(modewithdist[2] |> getpos, transformedphysicalmodesnext, 0.1, blockedcrystal) for modewithdist in transformedmodebydist[2]))

    courierseeds1RG = sum([reduce(+ ,FockMap(v,inspace = v |> getinspace |> orderedmodes |> setattr(:dumind => n) |> FockSpace,performpermute = false) for (n, (_ , v)) in enumerate(findlocalseeds(statecorrelations = couriercorrelationnext, regionfock = courierseedingfockRG, 
    spectrumextractpredicate = v -> true, symmetry = identitytransform(2)))) for courierseedingfockRG in courierseedingfocks1RG])
    courierseeds2RG = sum([reduce(+ ,FockMap(v,inspace = v |> getinspace |> orderedmodes |> setattr(:dumind => n) |> FockSpace,performpermute = false) for (n, (_ , v)) in enumerate(findlocalseeds(statecorrelations = couriercorrelationnext, regionfock = courierseedingfockRG, 
    spectrumextractpredicate = v -> true, symmetry = identitytransform(2)))) for courierseedingfockRG in courierseedingfocks2RG])

    courierseedsRG = courierseeds1RG + courierseeds2RG
    extendedcourierseedsRG = idmap(transformedlocalfocknext, transformedlocalfocknext)[:, courierseedsRG |> getoutspace] * courierseedsRG
    

    # Wannierization for couriers
    wannierizedcouriersRG = _localwannierization(transformedlocalmodesdict[:courier], extendedcourierseedsRG)
    
    return wannierizedfrozensRG, wannierizedcouriersRG
end

function locaclRGthirdnext(center_pt)::Tuple{FockMap,FockMap}
    center = Point(center_pt, scaledtriangular)
    # inspecting local region spectrum
    transformedphysicalmodessecond = couriercorrelationRGnext |> getoutspace |> orderedmodes
    
    trsasnformedlocalregionsecond, trsasnformedlocalfocksecond = localregioninspection(center,transformedphysicalmodessecond, 4, blockedcrystal)
    visualize(trsasnformedlocalregionsecond, title="l", visualspace=euclidean(RealSpace, 2))
    
    transformedlocaleigspec = couriercorrelationRGnext[trsasnformedlocalfocksecond, trsasnformedlocalfocksecond] |> eigspech 
    visualize(transformedlocaleigspec)
    
    transformedlocalmodesdictsecond = localmodesgrouping(couriercorrelationRGnext[trsasnformedlocalfocksecond, trsasnformedlocalfocksecond], 0.006)
    transformedmodebydistsecond = groupmodesbydist(region = trsasnformedlocalregionsecond,regionfock = trsasnformedlocalfocksecond,center = center)
    
    # finding seeds for local frozen (6 modes at the center)
    frozenseedingfocksRGsecond = (frozenseedingfock for (_,frozenseedingfock) in (localregioninspection(modewithdist[2] |> getpos, transformedphysicalmodessecond, 0.1, blockedcrystal) for modewithdist in transformedmodebydistsecond[2]))
    frozenseedsRGsecond = sum([reduce(+ ,FockMap(v,inspace = v |> getinspace |> orderedmodes |> setattr(:dumind => n) |> FockSpace,performpermute = false) for (n, (_ , v)) in enumerate(findlocalseeds(statecorrelations = couriercorrelationRGnext, regionfock = frozenseedingfockRG, 
    spectrumextractpredicate = v -> true, symmetry = identitytransform(2)))) for frozenseedingfockRG in frozenseedingfocksRGsecond])

    extendedfrozenseedsRGsecond = idmap(trsasnformedlocalfocksecond, trsasnformedlocalfocksecond)[:,frozenseedsRGsecond |> getoutspace] * frozenseedsRGsecond

    wannierizedfrozensRG = _localwannierization(transformedlocalmodesdictsecond[:frozen], extendedfrozenseedsRGsecond)

    # finding seeds for local courier (18 (12+6) modes at "boundary")
    courierseedingfocksRGsecond = (courierseedingfock for (_,courierseedingfock) in (localregioninspection(modewithdist[2] |> getpos, transformedphysicalmodessecond, 0.1, blockedcrystal) for modewithdist in transformedmodebydistsecond[1]))
    
    courierseedsRGsecond = sum([reduce(+ ,FockMap(v,inspace = v |> getinspace |> orderedmodes |> setattr(:dumind => n) |> FockSpace,performpermute = false) for (n, (_ , v)) in enumerate(findlocalseeds(statecorrelations = couriercorrelationRGnext, regionfock = courierseedingfockRG, 
    spectrumextractpredicate = v -> true, symmetry = identitytransform(2)))) for courierseedingfockRG in courierseedingfocksRGsecond])

    extendedcourierseedsRGsecond = idmap(trsasnformedlocalfocksecond, trsasnformedlocalfocksecond)[:, courierseedsRGsecond |> getoutspace] * courierseedsRGsecond
    

    # Wannierization for couriers
    wannierizedcouriersRG = _localwannierization(transformedlocalmodesdictsecond[:courier], extendedcourierseedsRGsecond)
    
    return wannierizedfrozensRG, wannierizedcouriersRG
end

wannier_fns_next = [locaclRGnext(center_pt) for center_pt in firstlayersnext]
localunitarylistnext = [wannierizedfrozen+wannierizedcourier for (wannierizedfrozen, wannierizedcourier) in wannier_fns_next]
localunitarynext = sum(localunitarylistnext)
wanniercourierlistnext = [wannierizedcourier for (_, wannierizedcourier) in wannier_fns_next]
wanniercouriersnext = sum(wanniercourierlistnext)

couriercorrelationnext = wanniercouriersnext'*couriercorrelationRGRG[localunitarynext |> getoutspace, localunitarynext |> getoutspace]*wanniercouriersnext


wannier_fnsRGnext = [locaclRGsecondnext(center_pt) for center_pt in secondlayersnext]
localunitarylistRGnext = [wannierizedfrozen+wannierizedcourier for (wannierizedfrozen, wannierizedcourier) in wannier_fnsRGnext]
localunitaryRGnext = sum(localunitarylistRGnext)
wanniercourierlistRGnext = [wannierizedcourier for (_, wannierizedcourier) in wannier_fnsRGnext]
wanniercouriersRGnext = sum(wanniercourierlistRGnext)

couriercorrelationRGnext = wanniercouriersRGnext'*couriercorrelationnext[localunitaryRGnext |> getoutspace, localunitaryRGnext |> getoutspace]*wanniercouriersRGnext

transformedlocalmodesdict = localmodesgrouping(couriercorrelationRGnext[trsasnformedlocalfocknext, trsasnformedlocalfocknext], 0.006)
transformedmodebydist = groupmodesbydist(region = trsasnformedlocalregionnext,regionfock = trsasnformedlocalfocknext,center = center)

wannier_fnsRGRGnext = [locaclRGthirdnext(center_pt) for center_pt in thirdlayersnext]
localunitarylistRGRGnext = [wannierizedfrozen+wannierizedcourier for (wannierizedfrozen, wannierizedcourier) in wannier_fnsRGRGnext]
localunitaryRGRGnext = sum(localunitarylistRGRGnext)
wanniercourierlistRGRGnext = [wannierizedcourier for (_, wannierizedcourier) in wannier_fnsRGRGnext]
wanniercouriersRGRGnext = sum(wanniercourierlistRGRGnext)

couriercorrelationRGRGnext = wanniercouriersRGRGnext'*couriercorrelationRGnext[localunitaryRGRGnext |> getoutspace, localunitaryRGRGnext |> getoutspace]*wanniercouriersRGRGnext

# inspecting local region spectrum
center_pt = [0,0]
center = Point(center_pt, scaledtriangular)
transformedphysicalmodesnext = couriercorrelationRGRGnext |> getoutspace |> orderedmodes

trsasnformedlocalregionnext, trsasnformedlocalfocknext = localregioninspection(center,transformedphysicalmodesnext, 4, blockedcrystal)
visualize(trsasnformedlocalregionnext, title="l", visualspace=euclidean(RealSpace, 2))

localeigspecnext = couriercorrelationRGRGnext[trsasnformedlocalfocknext, trsasnformedlocalfocknext] |> eigspech 
visualize(localeigspecnext)
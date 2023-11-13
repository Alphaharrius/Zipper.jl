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
        crystalspstates::Dict{Momentum, FockMap} = crystalisometries(localisometry=spstate, crystalfock=statecorrelations.outspace)
        crystalspstate::FockMap = directsum(v for (_, v) in crystalspstates)
        pseudoidentity::FockMap = (crystalspstate' * crystalspstate)
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
    seedingfock::FockSpace = FockSpace{Region}(seedingmodes)
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
        return inmode |> setattr(:offset => offset) |> setattr(:pos => basis) |> setattr(:ind => ind)
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
crystal = Crystal(unitcell, [48, 48])
reciprocalhashcalibration(crystal.sizes)

modes::Subset{Mode} = quantize(:pos, unitcell, 1)
m0, m1 = members(modes)

tₙ = ComplexF64(-1)
tₕ = 0.1im

nearestneighbor = [
    (m0, m1) => tₙ,
    (m0, setattr(m1, :offset => Point([-1, 0], triangular))) => tₙ,
    (m0, setattr(m1, :offset => Point([0, 1], triangular))) => tₙ]

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

C |> crystalspectrum |> visualize

correlations = C

crystalfock = correlations.outspace

scale = Scale([2 0; 0 2], crystalfock |> getcrystal |> getspace)
blockresult = blocking(:action => scale, :correlations => correlations, :crystal => crystalfock |> getcrystal)

blockedcrystal::Crystal = blockresult[:crystal]
blockedcorrelations::FockMap = blockresult[:correlations]

crystalpoints::Subset{Offset} = latticepoints(blockedcrystal)
blockedmodes::Subset{Mode} = quantize(:pos, blockedcrystal.unitcell, 1)
physicalmodes::Subset{Mode} = spanoffset(blockedmodes, crystalpoints)
scaledtriangular = scale*triangular 

center_pt = [-2,0]
center = Point(center_pt, scaledtriangular)

localregion,localfock = localregioninspection(center , physicalmodes, 2,blockedcrystal)
visualize(localregion)

localcorrelation = regioncorrelations(blockedcorrelations,localfock)
localeigspec = localcorrelation |> eigspech 
visualize(localeigspec)

modebydist = groupmodesbydist(region = localregion,regionfock = localfock,center = center)

# courierseedingfocks = (courierseedingfock for (courierseedingregion,courierseedingfock) in (localregioninspection(modewithdist[2] |> getpos, physicalmodes, 0.1, blockedcrystal) for modewithdist in modebydist[1]))


# courierseeds = [reduce(+ ,FockMap(v,inspace = v |> getinspace |> orderedmodes |> setattr(:dumind => n) |> FockSpace,performpermute = false) for (n, (_ , v)) in enumerate(_findlocalspstates(statecorrelations = blockedcorrelations, regionfock = courierseedingfock, 
#     spectrumextractpredicate = v -> true, symmetry = identitytransform(2)))) for courierseedingfock in courierseedingfocks]

# sum(courierseeds) |> getoutspace

# for modewithdist in modebydist[1]
#     localregioninspection(modewithdist[2] |> getpos, physicalmodes, 0.1, blockedcrystal)
# end


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

center_pt1 = [-4,-3]
wannierizedfrozens1, wannierizedcouriers1 = locaclRG(center_pt1)

center_pt2 = [-3,-1]
wannierizedfrozens2, wannierizedcouriers2 = locaclRG(center_pt2)

center_pt3 = [-2,1]
wannierizedfrozens3, wannierizedcouriers3 = locaclRG(center_pt3)

center_pt4 = [-1,3]
wannierizedfrozens4, wannierizedcouriers4 = locaclRG(center_pt4)

center_pt5 = [-3,-4]
wannierizedfrozens5, wannierizedcouriers5 = locaclRG(center_pt5)

center_pt6 = [-2,-2]
wannierizedfrozens6, wannierizedcouriers6 = locaclRG(center_pt6)

center_pt7 = [-1,0]
wannierizedfrozens7, wannierizedcouriers7 = locaclRG(center_pt7)

center_pt8 = [0,2]
wannierizedfrozens8, wannierizedcouriers8 = locaclRG(center_pt8)

center_pt9 = [1,4]
wannierizedfrozens9, wannierizedcouriers9 = locaclRG(center_pt9)

center_pt10 = [-1,-3]
wannierizedfrozens10, wannierizedcouriers10 = locaclRG(center_pt10)

center_pt11 = [0,-1]
wannierizedfrozens11, wannierizedcouriers11 = locaclRG(center_pt11)

center_pt12 = [1,1]
wannierizedfrozens12, wannierizedcouriers12 = locaclRG(center_pt12)

center_pt13 = [2,3]
wannierizedfrozens13, wannierizedcouriers13 = locaclRG(center_pt13)

center_pt14 = [1,-2]
wannierizedfrozens14, wannierizedcouriers14 = locaclRG(center_pt14)

center_pt15 = [2,0]
wannierizedfrozens15, wannierizedcouriers15 = locaclRG(center_pt15)

center_pt16 = [3,2]
wannierizedfrozens16, wannierizedcouriers16 = locaclRG(center_pt16)

center_pt17 = [3,-1]
wannierizedfrozens17, wannierizedcouriers17 = locaclRG(center_pt17)

center_pt18 = [4,1]
wannierizedfrozens18, wannierizedcouriers18 = locaclRG(center_pt18)



localunitary1 = wannierizedfrozens1 + wannierizedcouriers1
localunitary2 = wannierizedfrozens2 + wannierizedcouriers2
localunitary3 = wannierizedfrozens3 + wannierizedcouriers3
localunitary4 = wannierizedfrozens4 + wannierizedcouriers4
localunitary5 = wannierizedfrozens5 + wannierizedcouriers5
localunitary6 = wannierizedfrozens6 + wannierizedcouriers6
localunitary7 = wannierizedfrozens7 + wannierizedcouriers7
localunitary8 = wannierizedfrozens8 + wannierizedcouriers8
localunitary9 = wannierizedfrozens9 + wannierizedcouriers9
localunitary10 = wannierizedfrozens10 + wannierizedcouriers10
localunitary11 = wannierizedfrozens11 + wannierizedcouriers11
localunitary12 = wannierizedfrozens12 + wannierizedcouriers12
localunitary13 = wannierizedfrozens13 + wannierizedcouriers13
localunitary14 = wannierizedfrozens14 + wannierizedcouriers14
localunitary15 = wannierizedfrozens15 + wannierizedcouriers15
localunitary16 = wannierizedfrozens16 + wannierizedcouriers16
localunitary17 = wannierizedfrozens17 + wannierizedcouriers17
localunitary18 = wannierizedfrozens18 + wannierizedcouriers18

localunitary = localunitary1 + localunitary2 + localunitary3 + localunitary4 + localunitary5 + localunitary6 + localunitary7 + localunitary8 + localunitary9 + localunitary10 + localunitary11 + localunitary12 + localunitary13 + localunitary14 + localunitary15 + localunitary16 + localunitary17 + localunitary18

wannierizedcouriers = wannierizedcouriers1 + wannierizedcouriers2 + wannierizedcouriers3 + wannierizedcouriers4 + wannierizedcouriers5 + wannierizedcouriers6 + wannierizedcouriers7 + wannierizedcouriers8 + wannierizedcouriers9 + wannierizedcouriers10 + wannierizedcouriers11 + wannierizedcouriers12 + wannierizedcouriers13 + wannierizedcouriers14 + wannierizedcouriers15 + wannierizedcouriers16 + wannierizedcouriers17 + wannierizedcouriers18

couriercorrelation = wannierizedcouriers'*regioncorrelations(blockedcorrelations,localunitary |> getoutspace)*wannierizedcouriers

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

center_pt = [0,0]
center = Point(center_pt, scaledtriangular)
# inspecting local region spectrum
transformedphysicalmodes = couriercorrelation |> getoutspace |> orderedmodes
    
trsasnformedlocalregion, trsasnformedlocalfock = localregioninspection(center ,transformedphysicalmodes, 2, blockedcrystal)
visualize(trsasnformedlocalregion, title="l", visualspace=euclidean(RealSpace, 2))
    
transformedlocaleigspec = couriercorrelation[trsasnformedlocalfock, trsasnformedlocalfock] |> eigspech 
visualize(transformedlocaleigspec)

transformedmodebydist = groupmodesbydist(region = trsasnformedlocalregion,regionfock = trsasnformedlocalfock,center = center)

center_ptA = [-3,-2]
wannierizedfrozensRGA, wannierizedcouriersRGA = locaclRGsecond(center_ptA)::Tuple{FockMap,FockMap}

center_ptB = [-2,0]
wannierizedfrozensRGB, wannierizedcouriersRGB = locaclRGsecond(center_ptB)::Tuple{FockMap,FockMap}

center_ptC = [-1,2]
wannierizedfrozensRGC, wannierizedcouriersRGC = locaclRGsecond(center_ptC)::Tuple{FockMap,FockMap}

center_ptD = [-2,-3]
wannierizedfrozensRGD, wannierizedcouriersRGD = locaclRGsecond(center_ptD)::Tuple{FockMap,FockMap}

center_ptE = [-1,-1]
wannierizedfrozensRGE, wannierizedcouriersRGE = locaclRGsecond(center_ptE)::Tuple{FockMap,FockMap}

center_ptF = [0,1]
wannierizedfrozensRGF, wannierizedcouriersRGF = locaclRGsecond(center_ptF)::Tuple{FockMap,FockMap}

center_ptG = [1,3]
wannierizedfrozensRGG, wannierizedcouriersRGG = locaclRGsecond(center_ptG)::Tuple{FockMap,FockMap}

center_ptH = [0,-2]
wannierizedfrozensRGH, wannierizedcouriersRGH = locaclRGsecond(center_ptH)::Tuple{FockMap,FockMap}

center_ptI = [1,0]
wannierizedfrozensRGI, wannierizedcouriersRGI = locaclRGsecond(center_ptI)::Tuple{FockMap,FockMap}

center_ptJ = [2,2]
wannierizedfrozensRGJ, wannierizedcouriersRGJ = locaclRGsecond(center_ptJ)::Tuple{FockMap,FockMap}

center_ptK = [2,-1]
wannierizedfrozensRGK, wannierizedcouriersRGK = locaclRGsecond(center_ptK)::Tuple{FockMap,FockMap}

center_ptL = [3,1]
wannierizedfrozensRGL, wannierizedcouriersRGL = locaclRGsecond(center_ptL)::Tuple{FockMap,FockMap}

wannierizedcouriersRGA[:,3] |> Zipper.columnspec |> visualize


localunitaryRGsecondA = wannierizedfrozensRGA + wannierizedcouriersRGA
localunitaryRGsecondB = wannierizedfrozensRGB + wannierizedcouriersRGB
localunitaryRGsecondC = wannierizedfrozensRGC + wannierizedcouriersRGC
localunitaryRGsecondD = wannierizedfrozensRGD + wannierizedcouriersRGD
localunitaryRGsecondE = wannierizedfrozensRGE + wannierizedcouriersRGE
localunitaryRGsecondF = wannierizedfrozensRGF + wannierizedcouriersRGF
localunitaryRGsecondG = wannierizedfrozensRGG + wannierizedcouriersRGG
localunitaryRGsecondH = wannierizedfrozensRGH + wannierizedcouriersRGH
localunitaryRGsecondI = wannierizedfrozensRGI + wannierizedcouriersRGI
localunitaryRGsecondJ = wannierizedfrozensRGJ + wannierizedcouriersRGJ
localunitaryRGsecondK = wannierizedfrozensRGK + wannierizedcouriersRGK
localunitaryRGsecondL = wannierizedfrozensRGL + wannierizedcouriersRGL

localunitaryRG = localunitaryRGsecondA + localunitaryRGsecondB + localunitaryRGsecondC + localunitaryRGsecondD + localunitaryRGsecondE + localunitaryRGsecondF + localunitaryRGsecondG + localunitaryRGsecondH + localunitaryRGsecondI + localunitaryRGsecondJ + localunitaryRGsecondK + localunitaryRGsecondL

wannierizedcouriersRG = wannierizedcouriersRGA + wannierizedcouriersRGB + wannierizedcouriersRGC + wannierizedcouriersRGD + wannierizedcouriersRGE + wannierizedcouriersRGF + wannierizedcouriersRGG + wannierizedcouriersRGH + wannierizedcouriersRGI + wannierizedcouriersRGJ + wannierizedcouriersRGK + wannierizedcouriersRGL

couriercorrelationRG = wannierizedcouriersRG'*couriercorrelation[localunitaryRG |> getoutspace, localunitaryRG |> getoutspace]*wannierizedcouriersRG

# inspecting local region spectrum
center_pt = [-2,-1]
center = Point(center_pt, scaledtriangular)
transformedphysicalmodessecond = couriercorrelationRG |> getoutspace |> orderedmodes

trsasnformedlocalregionsecond, trsasnformedlocalfocksecond = localregioninspection(center,transformedphysicalmodessecond, 2, blockedcrystal)
visualize(trsasnformedlocalregionsecond, title="l", visualspace=euclidean(RealSpace, 2))

transformedlocaleigspec = couriercorrelationRG[trsasnformedlocalfocksecond, trsasnformedlocalfocksecond] |> eigspech 
visualize(transformedlocaleigspec)

transformedlocalmodesdictsecond = localmodesgrouping(couriercorrelationRG[trsasnformedlocalfocksecond, trsasnformedlocalfocksecond], 0.002)

transformedmodebydistsecond = groupmodesbydist(region = trsasnformedlocalregionsecond,regionfock = trsasnformedlocalfocksecond,center = center)

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

center_pta = [-2,-1]
wannierizedfrozensRGa, wannierizedcouriersRGa = locaclRGthird(center_pta)::Tuple{FockMap,FockMap}

center_ptb = [-1,1]
wannierizedfrozensRGb, wannierizedcouriersRGb = locaclRGthird(center_ptb)::Tuple{FockMap,FockMap}

center_ptc = [-1,-2]
wannierizedfrozensRGc, wannierizedcouriersRGc = locaclRGthird(center_ptc)::Tuple{FockMap,FockMap}

center_ptd = [0,0]
wannierizedfrozensRGd, wannierizedcouriersRGd = locaclRGthird(center_ptd)::Tuple{FockMap,FockMap}

center_pte = [1,2]
wannierizedfrozensRGe, wannierizedcouriersRGe = locaclRGthird(center_pte)::Tuple{FockMap,FockMap}

center_ptf = [1,-1]
wannierizedfrozensRGf, wannierizedcouriersRGf = locaclRGthird(center_ptf)::Tuple{FockMap,FockMap}

center_ptg = [2,1]
wannierizedfrozensRGg, wannierizedcouriersRGg = locaclRGthird(center_ptg)::Tuple{FockMap,FockMap}

localunitaryRGseconda = wannierizedfrozensRGa + wannierizedcouriersRGa
localunitaryRGsecondb = wannierizedfrozensRGb + wannierizedcouriersRGb
localunitaryRGsecondc = wannierizedfrozensRGc + wannierizedcouriersRGc
localunitaryRGsecondd = wannierizedfrozensRGd + wannierizedcouriersRGd
localunitaryRGseconde = wannierizedfrozensRGe + wannierizedcouriersRGe
localunitaryRGsecondf = wannierizedfrozensRGf + wannierizedcouriersRGf
localunitaryRGsecondg = wannierizedfrozensRGg + wannierizedcouriersRGg

localunitaryRGsecond = localunitaryRGseconda + localunitaryRGsecondb + localunitaryRGsecondc + localunitaryRGsecondd + localunitaryRGseconde + localunitaryRGsecondf + localunitaryRGsecondg 

wannierizedcouriersRGsecond = wannierizedcouriersRGa + wannierizedcouriersRGb + wannierizedcouriersRGc + wannierizedcouriersRGd + wannierizedcouriersRGe + wannierizedcouriersRGf + wannierizedcouriersRGg 

couriercorrelationRGsecond = wannierizedcouriersRGsecond'*couriercorrelationRG[localunitaryRGsecond |> getoutspace, localunitaryRGsecond |> getoutspace]*wannierizedcouriersRGsecond

# inspecting local region spectrum
center_pt = [0,0]
center = Point(center_pt, scaledtriangular)
transformedphysicalmodesnext = couriercorrelationRGsecond |> getoutspace |> orderedmodes

trsasnformedlocalregionnext, trsasnformedlocalfocknext = localregioninspection(center,transformedphysicalmodesnext, 4, blockedcrystal)
visualize(trsasnformedlocalregionnext, title="l", visualspace=euclidean(RealSpace, 2))

transformedlocaleigspec = couriercorrelationRGsecond[trsasnformedlocalfocknext, trsasnformedlocalfocknext] |> eigspech 
visualize(transformedlocaleigspec)


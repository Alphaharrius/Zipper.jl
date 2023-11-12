using Plotly, SmithNormalForm, LinearAlgebra, OrderedCollections, SparseArrays, Combinatorics
using Zipper

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
    symmetricspstates = (state * symmetricmap(symmetry, state) for state in orthogonalspstates)
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
    symmetricspstates = (state * symmetricmap(symmetry, state) for state in orthogonalspstates)
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
    symmetricspstates = (state * symmetricmap(symmetry, state) for state in orthogonalspstates)
    spstates = (state * spatialmap(state)' for state in symmetricspstates)

    return (state |> getinspace |> dimension => state for state in spstates)
end

function circularfilter(mode::Mode, center::Offset, radius::Real = 1.5)::Bool
    return norm(((mode-center) |> getpos |> euclidean))<=radius
end

function circularregionmodes(center::Offset, physicalmodes::Subset{Mode}, radius::Number, crystal:: Crystal)::Subset{Mode}
    currentspace::RealSpace = correlations.outspace |> getcrystal |> getspace |> orthogonalspace
    crystalsize = crystal |> size
    periodicnorm =  m -> lineartransform(currentspace, takeperiodic(m |> getpos,crystalsize)) |> norm
     
    function takeperiodic(p, crystalsize)
        return (p |> getspace)&[min(vectno,24-vectno) for (size,vectno) in zip(crystalsize, vec(p))]
    end
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

∈(vec, space::AffineSpace)::Point = Point(vec |> collect, space)

triangular = RealSpace([sqrt(3)/2 -1/2; 0. 1.]')
kspace = convert(MomentumSpace, triangular)

pa = [1/3, 2/3] ∈ triangular
pb = [2/3, 1/3] ∈ triangular
pc = (pa + pb) / 2
spatialsnappingcalibration((pa, pb, pc))

c6 = pointgrouptransform([cos(π/3) -sin(π/3); sin(π/3) cos(π/3)])
c3 = c6^2

unitcell = Subset(pa, pb)
crystal = Crystal(unitcell, [48, 48])
reciprocalhashcalibration(crystal.sizes)

modes::Subset{Mode} = quantize(:pos, unitcell, 1)
m0, m1 = members(modes)

t_a = ComplexF64(-0.3)
t_b = ComplexF64(-0.7)
t_n = ComplexF64(-0.5)
bonds::FockMap = bondmap([
    (m0, m0) => t_a,
    (m1, m1) => t_b,
    (m0, m1) => t_n,
    (m0, setattr(m1, :offset => Point([-1, 0], triangular))) => t_n,
    (m0, setattr(m1, :offset => Point([0, 1], triangular))) => t_n])

energyspectrum = computeenergyspectrum(bonds, crystal=crystal)
energyspectrum |> visualize

groundstates::CrystalSpectrum = groundstatespectrum(energyspectrum, perunitcellfillings=1)
groundstates |> visualize
groundstateprojector = groundstates |> crystalprojector

C = idmap(groundstateprojector.outspace) - groundstateprojector
correlationspectrum = C |> crystalspectrum
visualize(correlationspectrum, title="Physical Correlation")

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

localregion,localfock = localregioninspection(Point([0,0], scaledtriangular) , physicalmodes, 2, blockedcrystal)
visualize(localregion)


localcorrelation = regioncorrelations(blockedcorrelations,localfock)
localeigspec = localcorrelation |> eigspech 
visualize(localeigspec)

function locaclRG(center_pt)::Tuple{FockMap,FockMap}
    center = Point(center_pt, scaledtriangular)
    localregion,localfock = localregioninspection(center , physicalmodes, 2,blockedcrystal)
    visualize(localregion)

    localcorrelation = regioncorrelations(blockedcorrelations,localfock)
    localeigspec = localcorrelation |> eigspech 
    visualize(localeigspec)

    localmodesdict = localmodesgrouping(localcorrelation, 0.01)

    c6recenter = c6 |> recenter(Point(center_pt, scaledtriangular))
    c3recenter =c3 |> recenter(Point(center_pt, scaledtriangular))

    # finding seeds for local frozen (6 modes at the corners)
    frozenseedingcenterA::Offset = Point(center_pt+[2/3, 1/3] ,  (blockedmodes |> getspace))
    frozenseedingregionA,frozenseedingfockA = localregioninspection(frozenseedingcenterA, physicalmodes, 0.1, blockedcrystal)
    visualize(frozenseedingregionA, title="Frozen Seeding Region A", visualspace=euclidean(RealSpace, 2))
    frozenseedingfockB = FockSpace{Region}(m for m in c6recenter*frozenseedingfockA |> getoutspace)

    regioncorrelations(blockedcorrelations,frozenseedingfockB) |> eigspech |>visualize

    # for A 
    frozenseedA =  findlocalspstates(statecorrelations = blockedcorrelations, regionfock = frozenseedingfockA, 
        spectrumextractpredicate = v -> true, symmetry = identitytransform(2))[1]
    # for B 
    frozenseedB =  findlocalspstates(statecorrelations = blockedcorrelations, regionfock = frozenseedingfockB, 
        spectrumextractpredicate = v -> true, symmetry = identitytransform(2))[1]

    # Summing up all the symmetry related for corners
    frozenseedAs = frozenseedA + (c3recenter*frozenseedA.outspace)*frozenseedA*(c3recenter*frozenseedA.inspace)' + (c3recenter^2*frozenseedA.outspace)*frozenseedA*(c3recenter^2*frozenseedA.inspace)'
    frozenseedBs = frozenseedB + (c3recenter*frozenseedB.outspace)*frozenseedB*(c3recenter*frozenseedB.inspace)' + (c3recenter^2*frozenseedB.outspace)*frozenseedB*(c3recenter^2*frozenseedB.inspace)'

    # finding seeds for local frozen (6 modes at the center)
    frozenseedingcenter::Offset = Point(center_pt ,  (blockedmodes |> getspace))
    frozenseedingregion,frozenseedingfock = localregioninspection(frozenseedingcenter, physicalmodes, 1, blockedcrystal)
    visualize(frozenseedingregion, title="Frozen Seeding Region", visualspace=euclidean(RealSpace, 2))

    regioncorrelations(blockedcorrelations,frozenseedingfock) |> eigspech |>visualize

    frozenseedsorig = reduce(+ ,FockMap(v,inspace = v |> getinspace |> orderedmodes |> setattr(:dumind => n) |> FockSpace,performpermute = false) for (n, (_ , v)) in enumerate(_findlocalspstates(statecorrelations = blockedcorrelations, regionfock = frozenseedingfock, 
    spectrumextractpredicate = v -> true, symmetry = identitytransform(2))))

    # Summing all the frozen seed at center and corners to form a local frozen seeds
    frozenseeds = frozenseedsorig+frozenseedAs+frozenseedBs
    extendedfrozenseeds = idmap(localfock, localfock)[:,frozenseeds |> getoutspace] * frozenseeds

    wannierizedfrozens = _localwannierization(localmodesdict[:frozen], extendedfrozenseeds)

    # finding seeds for local courier (12 modes at the outermost boundary)
    courierseedingcenterA1::Offset = Point(center_pt+[2/3+1/6, 1/3-1/6] , (blockedmodes |> getspace))
    courierseedingregionA1,courierseedingfockA1 = localregioninspection(courierseedingcenterA1, physicalmodes, 0.1, blockedcrystal)
    visualize(courierseedingregionA1, title="Courier Seeding Region A1", visualspace=euclidean(RealSpace, 2))

    courierseedingcenterA2::Offset = Point(center_pt+[2/3+1/6, 1/3-1/6+1/2] , (blockedmodes |> getspace))
    courierseedingregionA2,courierseedingfockA2 = localregioninspection(courierseedingcenterA2, physicalmodes, 0.1, blockedcrystal)
    visualize(courierseedingregionA2, title="Courier Seeding Region A1", visualspace=euclidean(RealSpace, 2))

    courierseedingfockB1 = FockSpace{Region}(m for m in c6recenter*courierseedingfockA1 |> getoutspace)
    courierseedingfockB2 = FockSpace{Region}(m for m in c6recenter*courierseedingfockA2 |> getoutspace)

    # for A1
    courierseedA1 =  findlocalspstates(statecorrelations = blockedcorrelations, regionfock = courierseedingfockA1, 
        spectrumextractpredicate = v -> true, symmetry = identitytransform(2))[1]
    # for A2
    courierseedA2 =  findlocalspstates(statecorrelations = blockedcorrelations, regionfock = courierseedingfockA2, 
        spectrumextractpredicate = v -> true, symmetry = identitytransform(2))[1]

    # for B1
    courierseedB1 =  findlocalspstates(statecorrelations = blockedcorrelations, regionfock = courierseedingfockB1, 
        spectrumextractpredicate = v -> true, symmetry = identitytransform(2))[1]
    # for B2
    courierseedB2 =  findlocalspstates(statecorrelations = blockedcorrelations, regionfock = courierseedingfockB2, 
        spectrumextractpredicate = v -> true, symmetry = identitytransform(2))[1]

    # courierseedA = courierseedA1 + courierseedA2
    # courierseedB = courierseedB1 + courierseedB2

    # Summing up all the symmetry related for A (boundary)
    courierseedA1s = courierseedA1 + (c3recenter*courierseedA1.outspace)*courierseedA1*(c3recenter*courierseedA1.inspace)' + (c3recenter^2*courierseedA1.outspace)*courierseedA1*(c3recenter^2*courierseedA1.inspace)'
    courierseedA2s = courierseedA2 + (c3recenter*courierseedA2.outspace)*courierseedA2*(c3recenter*courierseedA2.inspace)' + (c3recenter^2*courierseedA2.outspace)*courierseedA2*(c3recenter^2*courierseedA2.inspace)'
    courierseedAs = courierseedA1s + courierseedA2s

    # Summing up all the symmetry related for B (boundary)
    courierseedB1s = courierseedB1 + (c3recenter*courierseedB1.outspace)*courierseedB1*(c3recenter*courierseedB1.inspace)' + (c3recenter^2*courierseedB1.outspace)*courierseedB1*(c3recenter^2*courierseedB1.inspace)'
    courierseedB2s = courierseedB2 + (c3recenter*courierseedB2.outspace)*courierseedB2*(c3recenter*courierseedB2.inspace)' + (c3recenter^2*courierseedB2.outspace)*courierseedB2*(c3recenter^2*courierseedB2.inspace)'
    courierseedBs = courierseedB1s + courierseedB2s

    courierseeds = courierseedAs + courierseedBs

    extendedcourierseeds = idmap(localfock, localfock)[:, courierseeds |> getoutspace] * courierseeds

    wannierizedcouriers = _localwannierization(localmodesdict[:courier], extendedcourierseeds)

    return wannierizedfrozens, wannierizedcouriers
end

center1 = [-2,-2]
wannierizedfrozens1, wannierizedcouriers1 = locaclRG(center1)

wannierizedcouriers1[:,5] |> Zipper.columnspec |> visualize

center2 = [-1,0]
wannierizedfrozens2, wannierizedcouriers2 = locaclRG(center2)

center3 = [0,2]
wannierizedfrozens3, wannierizedcouriers3 = locaclRG(center3)

center4 = [1,1]
wannierizedfrozens4, wannierizedcouriers4 = locaclRG(center4)

center5 = [2,0]
wannierizedfrozens5, wannierizedcouriers5 = locaclRG(center5)

center6 = [0,-1]
wannierizedfrozens6, wannierizedcouriers6 = locaclRG(center6)

localunitary1 = wannierizedfrozens1 + wannierizedcouriers1
localunitary2 = wannierizedfrozens2 + wannierizedcouriers2
localunitary3 = wannierizedfrozens3 + wannierizedcouriers3
localunitary4 = wannierizedfrozens4 + wannierizedcouriers4
localunitary5 = wannierizedfrozens5 + wannierizedcouriers5
localunitary6 = wannierizedfrozens6 + wannierizedcouriers6

# Getting real space form of the correlation matrix
ft = fourier(blockedcorrelations |> getoutspace, spanoffset(blockedcorrelations |> getinspace |> unitcellfock |> orderedmodes, blockedcrystal |> latticepoints)|> FockSpace)/sqrt(blockedcrystal |> vol)
blockedcorrelationsRS = ft'*blockedcorrelations*ft

# local "RGed" correlation matrix in real space
extendediso = ((((blockedcorrelationsRS|> getoutspace) - (localunitary1 |> getoutspace) - (localunitary2 |> getoutspace) - (localunitary3 |> getoutspace) - (localunitary4 |> getoutspace) - (localunitary5 |> getoutspace) - (localunitary6 |> getoutspace)) |> idmap)
                 + wannierizedcouriers1 +  wannierizedcouriers2 +  wannierizedcouriers3 +  wannierizedcouriers4 +  wannierizedcouriers5 +  wannierizedcouriers6)

transformedblockedcorrelationsRS = extendediso'*blockedcorrelationsRS*extendediso

# local RG step after the first RG
function locaclRGsecond(center_pt)::Tuple{FockMap,FockMap}
    center = Point(center_pt, scaledtriangular)
    # inspecting local region spectrum
    transformedphysicalmodes = transformedblockedcorrelationsRS |> getoutspace |> orderedmodes

    trsasnformedlocalregion, trsasnformedlocalfock = localregioninspection(center,transformedphysicalmodes, 2, blockedcrystal)
    visualize(trsasnformedlocalregion, title="l", visualspace=euclidean(RealSpace, 2))

    transformedlocaleigspec = transformedblockedcorrelationsRS[trsasnformedlocalfock, trsasnformedlocalfock] |> eigspech 
    visualize(transformedlocaleigspec)


    transformedlocalmodesdict = localmodesgrouping(transformedblockedcorrelationsRS[trsasnformedlocalfock, trsasnformedlocalfock], 0.02)

    # c6recenter = c6 |> recenter(scaledtriangular&center)
    # c3recenter = c3 |> recenter(scaledtriangular&center)

    # finding seeds for local frozen (6 modes at the centers)
    frozenseedingcenter::Offset = center
    frozenregionRG, frozenseedingfockRG = localregioninspection(frozenseedingcenter, transformedphysicalmodes, 1, blockedcrystal)
    visualize(frozenregionRG, visualspace=euclidean(RealSpace, 2))

    frozenRGspec = transformedblockedcorrelationsRS[frozenseedingfockRG, frozenseedingfockRG] |> eigspech 
    visualize(frozenRGspec)

    frozenseedsRG = reduce(+ ,FockMap(v,inspace = v |> getinspace |> orderedmodes |> setattr(:dumind => n) |> FockSpace,performpermute = false) for (n, (_ , v)) in enumerate(findlocalseeds(statecorrelations = transformedblockedcorrelationsRS, regionfock = frozenseedingfockRG, 
        spectrumextractpredicate = v -> true, symmetry = identitytransform(2))))

    extendedfrozenseedsRG = idmap(trsasnformedlocalfock, trsasnformedlocalfock)[:,frozenseedsRG |> getoutspace] * frozenseedsRG

    # Wannierization for frozen
    wannierizedfrozensRG = localwannierization(transformedlocalmodesdict[:frozen], extendedfrozenseedsRG)

    wannierizedfrozensRG[:,1] |> Zipper.columnspec |> visualize


    # finding seeds for local courier (6 modes at the corners)
    courierseedingcenterA1::Offset = Point(center_pt+[0 + 1/6, -3/4 + 1/12] ,  scaledtriangular)
    courierregionA1RG, courierseedingfockA1RG = localregioninspection(courierseedingcenterA1, transformedphysicalmodes, 0.1, blockedcrystal)
    visualize(courierregionA1RG, visualspace=euclidean(RealSpace, 2))

    courierseedingcenterB1::Offset = Point(center_pt+[0 - 1/6, -3/4 - 1/12] ,  scaledtriangular)
    courierregionB1RG, courierseedingfockB1RG = localregioninspection(courierseedingcenterB1, transformedphysicalmodes, 0.1, blockedcrystal)
    visualize(courierregionB1RG, visualspace=euclidean(RealSpace, 2))

    courierA1RGspec = transformedblockedcorrelationsRS[courierseedingfockA1RG, courierseedingfockA1RG] |> eigspech 
    visualize(courierA1RGspec)

    courierB1RGspec = transformedblockedcorrelationsRS[courierseedingfockB1RG, courierseedingfockB1RG] |> eigspech 
    visualize(courierB1RGspec)

    courierseedA1RG = reduce(+ ,FockMap(v,inspace = v |> getinspace |> orderedmodes |> setattr(:dumind => n) |> FockSpace,performpermute = false) for (n, (_ , v)) in enumerate(findlocalseeds(statecorrelations = transformedblockedcorrelationsRS, regionfock = courierseedingfockA1RG, 
    spectrumextractpredicate = v -> true, symmetry = identitytransform(2))))

    courierseedB1RG = reduce(+ ,FockMap(v,inspace = v |> getinspace |> orderedmodes |> setattr(:dumind => n) |> FockSpace,performpermute = false) for (n, (_ , v)) in enumerate(findlocalseeds(statecorrelations = transformedblockedcorrelationsRS, regionfock = courierseedingfockB1RG, 
    spectrumextractpredicate = v -> true, symmetry = identitytransform(2))))

    courierseedingcenterB2::Offset = Point(center_pt+[-3/4 + 1/12, 0 + 1/6] ,  scaledtriangular)
    courierregionB2RG, courierseedingfockB2RG = localregioninspection(courierseedingcenterB2, transformedphysicalmodes, 0.1, blockedcrystal)
    courierB2RGspec = transformedblockedcorrelationsRS[courierseedingfockB2RG, courierseedingfockB2RG] |> eigspech 
    visualize(courierB2RGspec)

    courierseedingcenterA2::Offset = Point(center_pt+[-3/4 - 1/12, 0 - 1/6] ,  scaledtriangular)
    courierregionA2RG, courierseedingfockA2RG = localregioninspection(courierseedingcenterA2, transformedphysicalmodes, 0.1, blockedcrystal)
    courierA2RGspec = transformedblockedcorrelationsRS[courierseedingfockA2RG, courierseedingfockA2RG] |> eigspech 
    visualize(courierA2RGspec)

    courierseedA2RG = reduce(+ ,FockMap(v,inspace = v |> getinspace |> orderedmodes |> setattr(:dumind => n) |> FockSpace,performpermute = false) for (n, (_ , v)) in enumerate(findlocalseeds(statecorrelations = transformedblockedcorrelationsRS, regionfock = courierseedingfockA2RG, 
    spectrumextractpredicate = v -> true, symmetry = identitytransform(2))))

    courierseedB2RG = reduce(+ ,FockMap(v,inspace = v |> getinspace |> orderedmodes |> setattr(:dumind => n) |> FockSpace,performpermute = false) for (n, (_ , v)) in enumerate(findlocalseeds(statecorrelations = transformedblockedcorrelationsRS, regionfock = courierseedingfockB2RG, 
    spectrumextractpredicate = v -> true, symmetry = identitytransform(2))))

    courierseedingcenterB3::Offset = Point(center_pt+[3/4 + 1/12, 3/4 - 1/12] ,  scaledtriangular)
    courierregionB3RG, courierseedingfockB3RG = localregioninspection(courierseedingcenterB3, transformedphysicalmodes, 0.1, blockedcrystal)
    courierB3RGspec = transformedblockedcorrelationsRS[courierseedingfockB3RG, courierseedingfockB3RG] |> eigspech 
    visualize(courierB3RGspec)

    courierseedingcenterA3::Offset = Point(center_pt+[3/4 - 1/12, 3/4 + 1/12] ,  scaledtriangular)
    courierregionA3RG, courierseedingfockA3RG = localregioninspection(courierseedingcenterA3, transformedphysicalmodes, 0.1, blockedcrystal)
    courierA3RGspec = transformedblockedcorrelationsRS[courierseedingfockA3RG, courierseedingfockA3RG] |> eigspech 
    visualize(courierA3RGspec)

    courierseedA3RG = reduce(+ ,FockMap(v,inspace = v |> getinspace |> orderedmodes |> setattr(:dumind => n) |> FockSpace,performpermute = false) for (n, (_ , v)) in enumerate(findlocalseeds(statecorrelations = transformedblockedcorrelationsRS, regionfock = courierseedingfockA3RG, 
    spectrumextractpredicate = v -> true, symmetry = identitytransform(2))))

    courierseedB3RG = reduce(+ ,FockMap(v,inspace = v |> getinspace |> orderedmodes |> setattr(:dumind => n) |> FockSpace,performpermute = false) for (n, (_ , v)) in enumerate(findlocalseeds(statecorrelations = transformedblockedcorrelationsRS, regionfock = courierseedingfockB3RG, 
    spectrumextractpredicate = v -> true, symmetry = identitytransform(2))))

    courierseedRG = courierseedA1RG + courierseedB1RG + courierseedA2RG + courierseedB2RG + courierseedA3RG + courierseedB3RG

    extendedcourierseedsRG = idmap(trsasnformedlocalfock, trsasnformedlocalfock)[:,courierseedRG |> getoutspace] * courierseedRG

    # Wannierization for couriers
    wannierizedcouriersRG = _localwannierization(transformedlocalmodesdict[:courier], extendedcourierseedsRG)


    return wannierizedfrozensRG, wannierizedcouriersRG
end

transformedphysicalmodes = transformedblockedcorrelationsRS |> getoutspace |> orderedmodes

trsasnformedlocalregion, trsasnformedlocalfock = localregioninspection(Point(center_pt ,  scaledtriangular),transformedphysicalmodes, 2, blockedcrystal)
visualize(trsasnformedlocalregion, title="l", visualspace=euclidean(RealSpace, 2))

transformedlocaleigspec = transformedblockedcorrelationsRS[trsasnformedlocalfock, trsasnformedlocalfock] |> eigspech 

visualize(transformedlocaleigspec)

centerA = [1,0]
wannierizedfrozensRGA, wannierizedcouriersRGA = locaclRGsecond(centerA)::Tuple{FockMap,FockMap}

wannierizedcouriersRGA[:,3] |> Zipper.columnspec |> visualize

centerB = [0,1]
wannierizedfrozensRGB, wannierizedcouriersRGB = locaclRGsecond(centerB)::Tuple{FockMap,FockMap}

centerC = [-1,-1]
wannierizedfrozensRGC, wannierizedcouriersRGC = locaclRGsecond(centerC)::Tuple{FockMap,FockMap}


localunitaryRGsecondA = wannierizedfrozensRGA + wannierizedcouriersRGA
localunitaryRGsecondB = wannierizedfrozensRGB + wannierizedcouriersRGB
localunitaryRGsecondC = wannierizedfrozensRGC + wannierizedcouriersRGC

extendedisosecond = ((((transformedblockedcorrelationsRS|> getoutspace) - (localunitaryRGsecondA  |> getoutspace) - (localunitaryRGsecondB  |> getoutspace) - (localunitaryRGsecondC  |> getoutspace)) |> idmap)
                 + wannierizedcouriersRGA + wannierizedcouriersRGB + wannierizedcouriersRGC)

transformedblockedcorrelationsRSsecond = extendedisosecond'*transformedblockedcorrelationsRS*extendedisosecond


# inspecting local region spectrum
transformedphysicalmodessecond = transformedblockedcorrelationsRSsecond |> getoutspace |> orderedmodes

trsasnformedlocalregionsecond, trsasnformedlocalfocksecond = localregioninspection(Point([0,0], scaledtriangular),transformedphysicalmodessecond, 2, blockedcrystal)
visualize(trsasnformedlocalregionsecond, title="l", visualspace=euclidean(RealSpace, 2))

transformedlocaleigspec = transformedblockedcorrelationsRSsecond[trsasnformedlocalfocksecond, trsasnformedlocalfocksecond] |> eigspech 
visualize(transformedlocaleigspec)
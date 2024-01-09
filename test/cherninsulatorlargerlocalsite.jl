using Plotly, SmithNormalForm, LinearAlgebra, OrderedCollections, SparseArrays, Combinatorics
using Revise
using Zipper
findlocalspstates

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
        return inmode |> setattr(:offset => offset) |> setattr(:b => basis) |> setattr(:ind => ind)
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

unitcell = Subset(pa, pb)
crystal = Crystal(unitcell, [48, 48])
reciprocalhashcalibration(crystal.sizes)

modes::Subset{Mode} = quantize(:b, unitcell, 1)
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
blockedmodes::Subset{Mode} = quantize(:b, blockedcrystal.unitcell, 1)
physicalmodes::Subset{Mode} = spanoffset(blockedmodes, crystalpoints)
scaledtriangular = scale*triangular 

center_pt = [-3,-3]
center = Point(center_pt, scaledtriangular)

localregion,localfock = localregioninspection(center , physicalmodes, 2.6,blockedcrystal)
visualize(localregion)

localcorrelation = regioncorrelations(blockedcorrelations,localfock)
localeigspec = localcorrelation |> eigspech 
visualize(localeigspec)

function locaclRG(center_pt)::Tuple{FockMap,FockMap}
    center = Point(center_pt, scaledtriangular)

    localregion,localfock = localregioninspection(center , physicalmodes, 2.6,blockedcrystal)
    visualize(localregion)

    localcorrelation = regioncorrelations(blockedcorrelations,localfock)
    localeigspec = localcorrelation |> eigspech 
    visualize(localeigspec)

    localmodesdict = localmodesgrouping(localcorrelation, 0.001)

    c6 = pointgrouptransform([cos(π/3) -sin(π/3); sin(π/3) cos(π/3)])
    c3 = c6^2

    c6recenter = c6 |> recenter(Point(center_pt, scaledtriangular))
    c3recenter = c3 |> recenter(Point(center_pt, scaledtriangular))

    # finding seeds for local frozen (6 modes at the center)
    frozenseedingcenter::Offset = center
    frozenseedingregion,frozenseedingfock = localregioninspection(frozenseedingcenter, physicalmodes, 2, blockedcrystal)
    visualize(frozenseedingregion, title="Frozen Seeding Region", visualspace=euclidean(RealSpace, 2))

    regioncorrelations(blockedcorrelations,frozenseedingfock) |> eigspech |>visualize
    frozenseedsorig = reduce(+ ,FockMap(v,inspace = v |> getinspace |> orderedmodes |> setattr(:dumind => n) |> FockSpace,performpermute = false) for (n, (_ , v)) in enumerate(_findlocalspstates(statecorrelations = blockedcorrelations, regionfock = frozenseedingfock, 
    spectrumextractpredicate = v -> true, symmetry = identitytransform(2))))

    extendedfrozenseeds = idmap(localfock, localfock)[:,frozenseedsorig |> getoutspace] * frozenseedsorig

    wannierizedfrozens = _localwannierization(localmodesdict[:frozen], extendedfrozenseeds)

    # finding seeds for local courier (6 modes at the corners)
    courierseedingcenterA::Offset = Point(center_pt+[1/6-1/3, 1+1/3-1/6] ,  (blockedmodes |> getspace))
    courierseedingregionA,courierseedingfockA = localregioninspection(courierseedingcenterA, physicalmodes, 0.5, blockedcrystal)
    visualize(courierseedingregionA, title="Courier Seeding Region A", visualspace=euclidean(RealSpace, 2))

    courierseedingcenterB::Offset = Point(center_pt+[1/6, 1+1/3] ,  (blockedmodes |> getspace))
    courierseedingregionB,courierseedingfockB = localregioninspection(courierseedingcenterB, physicalmodes, 0.1, blockedcrystal)
    visualize(courierseedingregionB, title="Courier Seeding Region B", visualspace=euclidean(RealSpace, 2))

    courierseedingcenterC::Offset = Point(center_pt+[1/3, 1+1/6] ,  (blockedmodes |> getspace))
    courierseedingregionC,courierseedingfockC = localregioninspection(courierseedingcenterC, physicalmodes, 0.1, blockedcrystal)
    visualize(courierseedingregionC, title="Courier Seeding Region C", visualspace=euclidean(RealSpace, 2))

    courierseedingcenterD::Offset = Point(center_pt+[1/3+1/3, 1+1/6+1/6] ,  (blockedmodes |> getspace))
    courierseedingregionD,courierseedingfockD = localregioninspection(courierseedingcenterD, physicalmodes, 0.1, blockedcrystal)
    visualize(courierseedingregionD, title="Courier Seeding Region D", visualspace=euclidean(RealSpace, 2))

    courierseedingcenterE::Offset = Point(center_pt+[1/3+1/3-1/6+1/3, 1+1/6+1/6-1/3+1/6] ,  (blockedmodes |> getspace))
    courierseedingregionE,courierseedingfockE = localregioninspection(courierseedingcenterE, physicalmodes, 0.1, blockedcrystal)
    visualize(courierseedingregionE, title="Courier Seeding Region E", visualspace=euclidean(RealSpace, 2))

    # for A 
    courierseedA =  findlocalspstates(statecorrelations = blockedcorrelations, regionfock = courierseedingfockA, 
        spectrumextractpredicate = v -> true, symmetry = identitytransform(2))[1]
    # for B 
    courierseedB =  findlocalspstates(statecorrelations = blockedcorrelations, regionfock = courierseedingfockB, 
        spectrumextractpredicate = v -> true, symmetry = identitytransform(2))[1]
    # for C
    courierseedC =  findlocalspstates(statecorrelations = blockedcorrelations, regionfock = courierseedingfockC, 
        spectrumextractpredicate = v -> true, symmetry = identitytransform(2))[1]
    # for D
    courierseedD =  findlocalspstates(statecorrelations = blockedcorrelations, regionfock = courierseedingfockD, 
        spectrumextractpredicate = v -> true, symmetry = identitytransform(2))[1]
    # for E 
    courierseedE =  findlocalspstates(statecorrelations = blockedcorrelations, regionfock = courierseedingfockE, 
        spectrumextractpredicate = v -> true, symmetry = identitytransform(2))[1]

    # Summing up all the symmetry related for corners
    courierseedAs = courierseedA + (c6recenter*courierseedA.outspace)*courierseedA*(c6recenter*courierseedA.inspace)' + (c6recenter^2*courierseedA.outspace)*courierseedA*(c6recenter^2*courierseedA.inspace)' + (c6recenter^3*courierseedA.outspace)*courierseedA*(c6recenter^3*courierseedA.inspace)' + (c6recenter^4*courierseedA.outspace)*courierseedA*(c6recenter^4*courierseedA.inspace)' + (c6recenter^5*courierseedA.outspace)*courierseedA*(c6recenter^5*courierseedA.inspace)'
    courierseedBs = courierseedB + (c6recenter*courierseedB.outspace)*courierseedB*(c6recenter*courierseedB.inspace)' + (c6recenter^2*courierseedB.outspace)*courierseedB*(c6recenter^2*courierseedB.inspace)' + (c6recenter^3*courierseedB.outspace)*courierseedB*(c6recenter^3*courierseedB.inspace)' + (c6recenter^4*courierseedB.outspace)*courierseedB*(c6recenter^4*courierseedB.inspace)' + (c6recenter^5*courierseedB.outspace)*courierseedB*(c6recenter^5*courierseedB.inspace)'
    courierseedCs = courierseedC + (c6recenter*courierseedC.outspace)*courierseedC*(c6recenter*courierseedC.inspace)' + (c6recenter^2*courierseedC.outspace)*courierseedC*(c6recenter^2*courierseedC.inspace)' + (c6recenter^3*courierseedC.outspace)*courierseedC*(c6recenter^3*courierseedC.inspace)' + (c6recenter^4*courierseedC.outspace)*courierseedC*(c6recenter^4*courierseedC.inspace)' + (c6recenter^5*courierseedC.outspace)*courierseedC*(c6recenter^5*courierseedC.inspace)'
    courierseedDs = courierseedD + (c6recenter*courierseedD.outspace)*courierseedD*(c6recenter*courierseedD.inspace)' + (c6recenter^2*courierseedD.outspace)*courierseedD*(c6recenter^2*courierseedD.inspace)' + (c6recenter^3*courierseedD.outspace)*courierseedD*(c6recenter^3*courierseedD.inspace)' + (c6recenter^4*courierseedD.outspace)*courierseedD*(c6recenter^4*courierseedD.inspace)' + (c6recenter^5*courierseedD.outspace)*courierseedD*(c6recenter^5*courierseedD.inspace)'
    courierseedEs = courierseedE + (c6recenter*courierseedE.outspace)*courierseedE*(c6recenter*courierseedE.inspace)' + (c6recenter^2*courierseedE.outspace)*courierseedE*(c6recenter^2*courierseedE.inspace)' + (c6recenter^3*courierseedE.outspace)*courierseedE*(c6recenter^3*courierseedE.inspace)' + (c6recenter^4*courierseedE.outspace)*courierseedE*(c6recenter^4*courierseedE.inspace)' + (c6recenter^5*courierseedE.outspace)*courierseedE*(c6recenter^5*courierseedE.inspace)'

    courierseeds = courierseedAs + courierseedBs + courierseedCs + courierseedDs + courierseedEs

    extendedcourierseeds = idmap(localfock, localfock)[:, courierseeds |> getoutspace] * courierseeds

    wannierizedcouriers = _localwannierization(localmodesdict[:courier], extendedcourierseeds)
    return wannierizedfrozens, wannierizedcouriers
end

center_pt1 = [-3,-3]
wannierizedfrozens1, wannierizedcouriers1 = locaclRG(center_pt1)

center_pt2 = [-1.5,0]
wannierizedfrozens2, wannierizedcouriers2 = locaclRG(center_pt2)

center_pt3 = [0,3]
wannierizedfrozens3, wannierizedcouriers3 = locaclRG(center_pt3)

center_pt4 = [1.5,1.5]
wannierizedfrozens4, wannierizedcouriers4 = locaclRG(center_pt4)

center_pt5 = [3,0]
wannierizedfrozens5, wannierizedcouriers5 = locaclRG(center_pt5)

center_pt6 = [0,-1.5]
wannierizedfrozens6, wannierizedcouriers6 = locaclRG(center_pt6)

localunitary1 = wannierizedfrozens1 + wannierizedcouriers1
localunitary2 = wannierizedfrozens2 + wannierizedcouriers2
localunitary3 = wannierizedfrozens3 + wannierizedcouriers3
localunitary4 = wannierizedfrozens4 + wannierizedcouriers4
localunitary5 = wannierizedfrozens5 + wannierizedcouriers5
localunitary6 = wannierizedfrozens6 + wannierizedcouriers6

localunitary = localunitary1 + localunitary2 + localunitary3 + localunitary4 + localunitary5 + localunitary6
wannierizedcouriers = wannierizedcouriers1 + wannierizedcouriers2 + wannierizedcouriers3 + wannierizedcouriers4 + wannierizedcouriers5 + wannierizedcouriers6

localunitary |> getoutspace
wannierizedcouriers'*regioncorrelations(blockedcorrelations,localunitary |> getoutspace)*wannierizedcouriers

# Getting real space form of the correlation matrix
ft = fourier(blockedcorrelations |> getoutspace, spanoffset(blockedcorrelations |> getinspace |> unitcellfock |> orderedmodes, blockedcrystal |> latticepoints)|> FockSpace)/sqrt(blockedcrystal |> vol)
blockedcorrelationsRS = ft'*blockedcorrelations*ft

# local "RGed" correlation matrix in real space
extendediso = ((((blockedcorrelationsRS|> getoutspace) - (localunitary1 |> getoutspace) - (localunitary2 |> getoutspace) - (localunitary3 |> getoutspace) - (localunitary4 |> getoutspace) - (localunitary5 |> getoutspace) - (localunitary6 |> getoutspace)) |> idmap)
                 + wannierizedcouriers1 +  wannierizedcouriers2 +  wannierizedcouriers3 +  wannierizedcouriers4 +  wannierizedcouriers5 +  wannierizedcouriers6)

# transformedblockedcorrelationsRS = extendediso'*blockedcorrelationsRS*extendediso
transformedblockedcorrelationsRS = wannierizedcouriers'*regioncorrelations(blockedcorrelations,localunitary |> getoutspace)*wannierizedcouriers

center_pt = [1.5,0]
center = Point(center_pt, scaledtriangular)
# inspecting local region spectrum
transformedphysicalmodes = transformedblockedcorrelationsRS |> getoutspace |> orderedmodes
    
trsasnformedlocalregion, trsasnformedlocalfock = localregioninspection(center ,transformedphysicalmodes, 2.6, blockedcrystal)
visualize(trsasnformedlocalregion, title="l", visualspace=euclidean(RealSpace, 2))
    
transformedlocaleigspec = transformedblockedcorrelationsRS[trsasnformedlocalfock, trsasnformedlocalfock] |> eigspech 
visualize(transformedlocaleigspec)

function locaclRGsecond(center_pt)::Tuple{FockMap,FockMap}
    center = Point(center_pt, scaledtriangular)
    # inspecting local region spectrum
    transformedphysicalmodes = transformedblockedcorrelationsRS |> getoutspace |> orderedmodes
        
    trsasnformedlocalregion, trsasnformedlocalfock = localregioninspection(center ,transformedphysicalmodes, 2.6, blockedcrystal)
    visualize(trsasnformedlocalregion, title="l", visualspace=euclidean(RealSpace, 2))
        
    transformedlocaleigspec = transformedblockedcorrelationsRS[trsasnformedlocalfock, trsasnformedlocalfock] |> eigspech 
    visualize(transformedlocaleigspec)
        
    transformedlocalmodesdict = localmodesgrouping(transformedblockedcorrelationsRS[trsasnformedlocalfock, trsasnformedlocalfock], 0.006)
        
        # c6recenter = c6 |> recenter(scaledtriangular&center)
        # c3recenter = c3 |> recenter(scaledtriangular&center)
        
    # finding seeds for local frozen (6 modes at the centers)
    frozenseedingcenter::Offset = Point(center_pt, scaledtriangular)
    frozenregionRG, frozenseedingfockRG = localregioninspection(frozenseedingcenter, transformedphysicalmodes, 2, blockedcrystal)
    visualize(frozenregionRG, visualspace=euclidean(RealSpace, 2))
        
    frozenRGspec = transformedblockedcorrelationsRS[frozenseedingfockRG, frozenseedingfockRG] |> eigspech 
    visualize(frozenRGspec)
        
    frozenseedsRG = reduce(+ ,FockMap(v,inspace = v |> getinspace |> orderedmodes |> setattr(:dumind => n) |> FockSpace,performpermute = false) for (n, (_ , v)) in enumerate(findlocalseeds(statecorrelations = transformedblockedcorrelationsRS, regionfock = frozenseedingfockRG, 
        spectrumextractpredicate = v -> true, symmetry = identitytransform(2))))
        
    extendedfrozenseedsRG = idmap(trsasnformedlocalfock, trsasnformedlocalfock)[:,frozenseedsRG |> getoutspace] * frozenseedsRG
        
    # Wannierization for frozen
    wannierizedfrozensRG = localwannierization(transformedlocalmodesdict[:frozen], extendedfrozenseedsRG)
        
    wannierizedfrozensRG[:,12] |> Zipper.columnspec |> visualize
        
        
    # finding seeds for local courier (6 modes closer to the center)
    courierseedingcenternA1::Offset = Point(center_pt+[1+1/6-1/3, 1+1/3-1/6] ,  scaledtriangular)
    courierseedingregionnA1RG,courierseedingfocknA1RG = localregioninspection(courierseedingcenternA1, transformedphysicalmodes, 0.1, blockedcrystal)
    visualize(courierseedingregionnA1RG, title="Courier Seeding Region A", visualspace=euclidean(RealSpace, 2))
        
    courierseedingcenternB1::Offset = Point(center_pt+[1+1/6, 1+1/3] ,  scaledtriangular)
    courierseedingregionnB1RG,courierseedingfocknB1RG = localregioninspection(courierseedingcenternB1, transformedphysicalmodes, 0.1, blockedcrystal)
    visualize(courierseedingregionnB1RG, title="Courier Seeding Region B", visualspace=euclidean(RealSpace, 2))

    courierseedingcenternC1::Offset = Point(center_pt+[1+1/3, 1+1/6] ,  scaledtriangular)
    courierseedingregionnC1RG,courierseedingfocknC1RG = localregioninspection(courierseedingcenternC1, transformedphysicalmodes, 0.1, blockedcrystal)
    visualize(courierseedingregionnC1RG, title="Courier Seeding Region C", visualspace=euclidean(RealSpace, 2))
        
    courierseedingcenternD1::Offset = Point(center_pt+[1-1/6+1/3, 1-1/3+1/6] ,  scaledtriangular)
    courierseedingregionnD1RG,courierseedingfocknD1RG = localregioninspection(courierseedingcenternD1, transformedphysicalmodes, 0.1, blockedcrystal)
    visualize(courierseedingregionnD1RG, title="Courier Seeding Region D", visualspace=euclidean(RealSpace, 2))
        
    courierseednA1 = reduce(+ ,FockMap(v,inspace = v |> getinspace |> orderedmodes |> setattr(:dumind => n) |> FockSpace,performpermute = false) for (n, (_ , v)) in enumerate(findlocalseeds(statecorrelations = transformedblockedcorrelationsRS, regionfock = courierseedingfocknA1RG, 
        spectrumextractpredicate = v -> true, symmetry = identitytransform(2))))
    courierseednB1 = reduce(+ ,FockMap(v,inspace = v |> getinspace |> orderedmodes |> setattr(:dumind => n) |> FockSpace,performpermute = false) for (n, (_ , v)) in enumerate(findlocalseeds(statecorrelations = transformedblockedcorrelationsRS, regionfock = courierseedingfocknB1RG, 
        spectrumextractpredicate = v -> true, symmetry = identitytransform(2))))
    courierseednC1 = reduce(+ ,FockMap(v,inspace = v |> getinspace |> orderedmodes |> setattr(:dumind => n) |> FockSpace,performpermute = false) for (n, (_ , v)) in enumerate(findlocalseeds(statecorrelations = transformedblockedcorrelationsRS, regionfock = courierseedingfocknC1RG, 
        spectrumextractpredicate = v -> true, symmetry = identitytransform(2))))
    courierseednD1 = reduce(+ ,FockMap(v,inspace = v |> getinspace |> orderedmodes |> setattr(:dumind => n) |> FockSpace,performpermute = false) for (n, (_ , v)) in enumerate(findlocalseeds(statecorrelations = transformedblockedcorrelationsRS, regionfock = courierseedingfocknD1RG, 
        spectrumextractpredicate = v -> true, symmetry = identitytransform(2))))

        
        
    courierseedingcenternA2::Offset = Point(center_pt+[0-1/3, -1-1/6] ,  scaledtriangular)
    courierseedingregionnA2RG,courierseedingfocknA2RG = localregioninspection(courierseedingcenternA2, transformedphysicalmodes, 0.1, blockedcrystal)
    visualize(courierseedingregionnA2RG, title="Courier Seeding Region A", visualspace=euclidean(RealSpace, 2))
        
    courierseedingcenternB2::Offset = Point(center_pt+[0-1/6, -1-1/3] ,  scaledtriangular)
    courierseedingregionnB2RG,courierseedingfocknB2RG = localregioninspection(courierseedingcenternB2, transformedphysicalmodes, 0.1, blockedcrystal)
    visualize(courierseedingregionnB2RG, title="Courier Seeding Region B", visualspace=euclidean(RealSpace, 2))

    courierseedingcenternC2::Offset = Point(center_pt+[0-1/6+1/3, -1-1/3+1/6],  scaledtriangular)
    courierseedingregionnC2RG,courierseedingfocknC2RG = localregioninspection(courierseedingcenternC2, transformedphysicalmodes, 0.1, blockedcrystal)
    visualize(courierseedingregionnC2RG, title="Courier Seeding Region C", visualspace=euclidean(RealSpace, 2))
        
    courierseedingcenternD2::Offset = Point(center_pt+[0+1/3, -1+1/6],  scaledtriangular)
    courierseedingregionnD2RG,courierseedingfocknD2RG = localregioninspection(courierseedingcenternD2, transformedphysicalmodes, 0.1, blockedcrystal)
    visualize(courierseedingregionnD2RG, title="Courier Seeding Region D", visualspace=euclidean(RealSpace, 2))
        
    courierseednA2 = reduce(+ ,FockMap(v,inspace = v |> getinspace |> orderedmodes |> setattr(:dumind => n) |> FockSpace,performpermute = false) for (n, (_ , v)) in enumerate(findlocalseeds(statecorrelations = transformedblockedcorrelationsRS, regionfock = courierseedingfocknA2RG, 
        spectrumextractpredicate = v -> true, symmetry = identitytransform(2))))
    courierseednB2 = reduce(+ ,FockMap(v,inspace = v |> getinspace |> orderedmodes |> setattr(:dumind => n) |> FockSpace,performpermute = false) for (n, (_ , v)) in enumerate(findlocalseeds(statecorrelations = transformedblockedcorrelationsRS, regionfock = courierseedingfocknB2RG, 
        spectrumextractpredicate = v -> true, symmetry = identitytransform(2))))
    courierseednC2 = reduce(+ ,FockMap(v,inspace = v |> getinspace |> orderedmodes |> setattr(:dumind => n) |> FockSpace,performpermute = false) for (n, (_ , v)) in enumerate(findlocalseeds(statecorrelations = transformedblockedcorrelationsRS, regionfock = courierseedingfocknC2RG, 
        spectrumextractpredicate = v -> true, symmetry = identitytransform(2))))
    courierseednD2 = reduce(+ ,FockMap(v,inspace = v |> getinspace |> orderedmodes |> setattr(:dumind => n) |> FockSpace,performpermute = false) for (n, (_ , v)) in enumerate(findlocalseeds(statecorrelations = transformedblockedcorrelationsRS, regionfock = courierseedingfocknD2RG, 
        spectrumextractpredicate = v -> true, symmetry = identitytransform(2))))

    courierseedingcenternA3::Offset = Point(center_pt+[-1+1/6, 0+1/3] ,  scaledtriangular)
    courierseedingregionnA3RG,courierseedingfocknA3RG = localregioninspection(courierseedingcenternA3, transformedphysicalmodes, 0.1, blockedcrystal)
    visualize(courierseedingregionnA3RG, title="Courier Seeding Region A", visualspace=euclidean(RealSpace, 2))
        
    courierseedingcenternB3::Offset = Point(center_pt+[-1+1/6-1/3, 0+1/3-1/6] ,  scaledtriangular)
    courierseedingregionnB3RG,courierseedingfocknB3RG = localregioninspection(courierseedingcenternB3, transformedphysicalmodes, 0.1, blockedcrystal)
    visualize(courierseedingregionnB3RG, title="Courier Seeding Region B", visualspace=euclidean(RealSpace, 2))

    courierseedingcenternC3::Offset = Point(center_pt+[-1-1/6, 0-1/3] ,  scaledtriangular)
    courierseedingregionnC3RG,courierseedingfocknC3RG = localregioninspection(courierseedingcenternC3, transformedphysicalmodes, 0.1, blockedcrystal)
    visualize(courierseedingregionnC3RG, title="Courier Seeding Region C", visualspace=euclidean(RealSpace, 2))
        
    courierseedingcenternD3::Offset = Point(center_pt+[-1-1/3, 0-1/6] ,  scaledtriangular)
    courierseedingregionnD3RG,courierseedingfocknD3RG = localregioninspection(courierseedingcenternD3, transformedphysicalmodes, 0.1, blockedcrystal)
    visualize(courierseedingregionnD3RG, title="Courier Seeding Region D", visualspace=euclidean(RealSpace, 2))
        
    courierseednA3 = reduce(+ ,FockMap(v,inspace = v |> getinspace |> orderedmodes |> setattr(:dumind => n) |> FockSpace,performpermute = false) for (n, (_ , v)) in enumerate(findlocalseeds(statecorrelations = transformedblockedcorrelationsRS, regionfock = courierseedingfocknA3RG, 
        spectrumextractpredicate = v -> true, symmetry = identitytransform(2))))
    courierseednB3 = reduce(+ ,FockMap(v,inspace = v |> getinspace |> orderedmodes |> setattr(:dumind => n) |> FockSpace,performpermute = false) for (n, (_ , v)) in enumerate(findlocalseeds(statecorrelations = transformedblockedcorrelationsRS, regionfock = courierseedingfocknB3RG, 
        spectrumextractpredicate = v -> true, symmetry = identitytransform(2))))
    courierseednC3 = reduce(+ ,FockMap(v,inspace = v |> getinspace |> orderedmodes |> setattr(:dumind => n) |> FockSpace,performpermute = false) for (n, (_ , v)) in enumerate(findlocalseeds(statecorrelations = transformedblockedcorrelationsRS, regionfock = courierseedingfocknC3RG, 
        spectrumextractpredicate = v -> true, symmetry = identitytransform(2))))
    courierseednD3 = reduce(+ ,FockMap(v,inspace = v |> getinspace |> orderedmodes |> setattr(:dumind => n) |> FockSpace,performpermute = false) for (n, (_ , v)) in enumerate(findlocalseeds(statecorrelations = transformedblockedcorrelationsRS, regionfock = courierseedingfocknD3RG, 
        spectrumextractpredicate = v -> true, symmetry = identitytransform(2))))
        
    courierseedRG = courierseednA1 + courierseednB1 + courierseednC1 + courierseednD1 + courierseednA2 + courierseednB2 + courierseednC2 + courierseednD2 + courierseednA3+ courierseednB3 + courierseednC3 + courierseednD3
        
    extendedcourierseedsRG = idmap(trsasnformedlocalfock, trsasnformedlocalfock)[:,courierseedRG |> getoutspace] * courierseedRG
        
    # Wannierization for couriers
    wannierizedcouriersRG = _localwannierization(transformedlocalmodesdict[:courier], extendedcourierseedsRG)
    return wannierizedfrozensRG, wannierizedcouriersRG
end
    

center_ptA = [1.5,0]
wannierizedfrozensRGA, wannierizedcouriersRGA = locaclRGsecond(center_ptA)::Tuple{FockMap,FockMap}

wannierizedcouriersRGA[:,10] |> Zipper.columnspec |> visualize

center_ptB = [0,1.5]
wannierizedfrozensRGB, wannierizedcouriersRGB = locaclRGsecond(center_ptB)::Tuple{FockMap,FockMap}

center_ptC = [-1.5,-1.5]
wannierizedfrozensRGC, wannierizedcouriersRGC = locaclRGsecond(center_ptC)::Tuple{FockMap,FockMap}


localunitaryRGsecondA = wannierizedfrozensRGA + wannierizedcouriersRGA
localunitaryRGsecondB = wannierizedfrozensRGB + wannierizedcouriersRGB
localunitaryRGsecondC = wannierizedfrozensRGC + wannierizedcouriersRGC

localunitaryRGsecond = localunitaryRGsecondA + localunitaryRGsecondB  + localunitaryRGsecondC
wannierizedcouriersRGsecond = wannierizedcouriersRGA + wannierizedcouriersRGB + wannierizedcouriersRGC

localunitaryRGsecond |> getoutspace
wannierizedcouriersRGsecond'*regioncorrelations(transformedblockedcorrelationsRS,localunitaryRGsecond |> getoutspace)*wannierizedcouriersRGsecond

extendedisosecond = ((((transformedblockedcorrelationsRS|> getoutspace) - (localunitaryRGsecondA  |> getoutspace) - (localunitaryRGsecondB  |> getoutspace) - (localunitaryRGsecondC  |> getoutspace)) |> idmap)
                 + wannierizedcouriersRGA + wannierizedcouriersRGB + wannierizedcouriersRGC)

# transformedblockedcorrelationsRSsecond = extendedisosecond'*transformedblockedcorrelationsRS*extendedisosecond
transformedblockedcorrelationsRSsecond = wannierizedcouriersRGsecond'*regioncorrelations(transformedblockedcorrelationsRS,localunitaryRGsecond |> getoutspace)*wannierizedcouriersRGsecond

# inspecting local region spectrum
center_pt = [0,0]
transformedphysicalmodessecond = transformedblockedcorrelationsRSsecond |> getoutspace |> orderedmodes

trsasnformedlocalregionsecond, trsasnformedlocalfocksecond = localregioninspection(Point([0,0], scaledtriangular),transformedphysicalmodessecond, 2, blockedcrystal)
visualize(trsasnformedlocalregionsecond, title="l", visualspace=euclidean(RealSpace, 2))

transformedlocaleigspec = transformedblockedcorrelationsRSsecond[trsasnformedlocalfocksecond, trsasnformedlocalfocksecond] |> eigspech 
visualize(transformedlocaleigspec)

transformedlocalmodesdictsecond = localmodesgrouping(transformedblockedcorrelationsRSsecond[trsasnformedlocalfocksecond, trsasnformedlocalfocksecond], 0.001)

# finding seeds for local frozen (6 modes at the centers)
frozenseedingcenter::Offset = Point(center_pt, scaledtriangular)
frozenregionRG, frozenseedingfockRG = localregioninspection(frozenseedingcenter, transformedphysicalmodessecond , 1, blockedcrystal)
visualize(frozenregionRG, visualspace=euclidean(RealSpace, 2))

frozenRGspec = transformedblockedcorrelationsRSsecond[frozenseedingfockRG, frozenseedingfockRG] |> eigspech 
visualize(frozenRGspec)

frozenseedsRG = reduce(+ ,FockMap(v,inspace = v |> getinspace |> orderedmodes |> setattr(:dumind => n) |> FockSpace,performpermute = false) for (n, (_ , v)) in enumerate(findlocalseeds(statecorrelations = transformedblockedcorrelationsRSsecond, regionfock = frozenseedingfockRG, 
    spectrumextractpredicate = v -> true, symmetry = identitytransform(2))))

extendedfrozenseedsRG = idmap(trsasnformedlocalfocksecond, trsasnformedlocalfocksecond)[:,frozenseedsRG |> getoutspace] * frozenseedsRG

# Wannierization for frozen
wannierizedfrozensRG = localwannierization(transformedlocalmodesdictsecond[:frozen], extendedfrozenseedsRG)

wannierizedfrozensRG[:,1] |> Zipper.columnspec |> visualize


# finding seeds for local courier (6 modes closer to the center)
courierseedingcenternA1::Offset = Point(center_pt+[2/3, 1/3] ,  scaledtriangular)
courierseedingregionnA1RG,courierseedingfocknA1RG = localregioninspection(courierseedingcenternA1, transformedphysicalmodessecond, 0.1, blockedcrystal)
# visualize(courierseedingregionnA1, title="Courier Seeding Region A", visualspace=euclidean(RealSpace, 2))

courierseedingcenternB1::Offset = Point(center_pt+[1/3, 2/3] ,  scaledtriangular)
courierseedingregionnB1RG,courierseedingfocknB1RG = localregioninspection(courierseedingcenternB1, transformedphysicalmodessecond, 0.1, blockedcrystal)
# visualize(courierseedingregionnB1, title="Courier Seeding Region B", visualspace=euclidean(RealSpace, 2))

# for A 
courierseednA1 = reduce(+ ,FockMap(v,inspace = v |> getinspace |> orderedmodes |> setattr(:dumind => n) |> FockSpace,performpermute = false) for (n, (_ , v)) in enumerate(findlocalseeds(statecorrelations = transformedblockedcorrelationsRSsecond, regionfock = courierseedingfocknA1RG, 
spectrumextractpredicate = v -> true, symmetry = identitytransform(2))))

# for B 
courierseednB1 =  reduce(+ ,FockMap(v,inspace = v |> getinspace |> orderedmodes |> setattr(:dumind => n) |> FockSpace,performpermute = false) for (n, (_ , v)) in enumerate(findlocalseeds(statecorrelations = transformedblockedcorrelationsRSsecond, regionfock = courierseedingfocknB1RG, 
spectrumextractpredicate = v -> true, symmetry = identitytransform(2))))


courierseedingcenternA2::Offset = Point(center_pt+[-1/3, -2/3] ,  scaledtriangular)
courierseedingregionnA2RG,courierseedingfocknA2RG = localregioninspection(courierseedingcenternA2, transformedphysicalmodessecond, 0.1, blockedcrystal)
visualize(courierseedingregionnA2RG, title="Courier Seeding Region A", visualspace=euclidean(RealSpace, 2))

courierseedingcenternB2::Offset = Point(center_pt+[1/3, -1/3] ,  scaledtriangular)
courierseedingregionnB2RG,courierseedingfocknB2RG = localregioninspection(courierseedingcenternB2, transformedphysicalmodessecond, 0.1, blockedcrystal)
visualize(courierseedingregionnB2RG, title="Courier Seeding Region B", visualspace=euclidean(RealSpace, 2))

# for A 
courierseednA2 =  reduce(+ ,FockMap(v,inspace = v |> getinspace |> orderedmodes |> setattr(:dumind => n) |> FockSpace,performpermute = false) for (n, (_ , v)) in enumerate(findlocalseeds(statecorrelations = transformedblockedcorrelationsRSsecond, regionfock = courierseedingfocknA2RG, 
spectrumextractpredicate = v -> true, symmetry = identitytransform(2))))
# for B 
courierseednB2 =  reduce(+ ,FockMap(v,inspace = v |> getinspace |> orderedmodes |> setattr(:dumind => n) |> FockSpace,performpermute = false) for (n, (_ , v)) in enumerate(findlocalseeds(statecorrelations = transformedblockedcorrelationsRSsecond, regionfock = courierseedingfocknB2RG, 
spectrumextractpredicate = v -> true, symmetry = identitytransform(2))))


courierseedingcenternA3::Offset = Point(center_pt+[-1/3, 1/3] ,  scaledtriangular)
courierseedingregionnA3RG,courierseedingfocknA3RG = localregioninspection(courierseedingcenternA3, transformedphysicalmodessecond, 0.1, blockedcrystal)
visualize(courierseedingregionnA2RG, title="Courier Seeding Region A", visualspace=euclidean(RealSpace, 2))

courierseedingcenternB3::Offset = Point(center_pt+[-2/3, -1/3] ,  scaledtriangular)
courierseedingregionnB3RG,courierseedingfocknB3RG = localregioninspection(courierseedingcenternB3, transformedphysicalmodessecond, 0.1, blockedcrystal)
visualize(courierseedingregionnB2RG, title="Courier Seeding Region B", visualspace=euclidean(RealSpace, 2))

# for A 
courierseednA3 =  reduce(+ ,FockMap(v,inspace = v |> getinspace |> orderedmodes |> setattr(:dumind => n) |> FockSpace,performpermute = false) for (n, (_ , v)) in enumerate(findlocalseeds(statecorrelations = transformedblockedcorrelationsRSsecond, regionfock = courierseedingfocknA3RG, 
spectrumextractpredicate = v -> true, symmetry = identitytransform(2))))
# for B 
courierseednB3 =  reduce(+ ,FockMap(v,inspace = v |> getinspace |> orderedmodes |> setattr(:dumind => n) |> FockSpace,performpermute = false) for (n, (_ , v)) in enumerate(findlocalseeds(statecorrelations = transformedblockedcorrelationsRSsecond, regionfock = courierseedingfocknB3RG, 
spectrumextractpredicate = v -> true, symmetry = identitytransform(2))))

courierseedRG = courierseednA1 + courierseednB1 + courierseednA2 + courierseednB2 + courierseednA3+ courierseednB3
    
extendedcourierseedsRG = idmap(trsasnformedlocalfocksecond, trsasnformedlocalfocksecond)[:,courierseedRG |> getoutspace] * courierseedRG
    
# Wannierization for couriers
wannierizedcouriersRG = _localwannierization(transformedlocalmodesdictsecond[:courier], extendedcourierseedsRG)

wannierizedcouriersRG[:,2] |> Zipper.columnspec |> visualize

localunitaryRGthird = wannierizedcouriersRG + wannierizedfrozensRG

extendedisothird = ((((transformedblockedcorrelationsRSsecond|> getoutspace) - (localunitaryRGthird  |> getoutspace)) |> idmap)
                 + wannierizedcouriersRG)

transformedblockedcorrelationsRSthird = extendedisothird'*transformedblockedcorrelationsRSsecond*extendedisothird

center_pt = [0,0]
transformedphysicalmodesthird = transformedblockedcorrelationsRSthird |> getoutspace |> orderedmodes

trsasnformedlocalregionthird, trsasnformedlocalfockthird = localregioninspection(Point([0,0], scaledtriangular),transformedphysicalmodesthird, 2, blockedcrystal)
visualize(trsasnformedlocalregionthird, title="l", visualspace=euclidean(RealSpace, 2))

transformedlocaleigspec = transformedblockedcorrelationsRSthird[trsasnformedlocalfockthird, trsasnformedlocalfockthird] |> eigspech 
visualize(transformedlocaleigspec)
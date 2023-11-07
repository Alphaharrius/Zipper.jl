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

function circularregionmodes(center::Offset, physicalmodes::Subset{Mode}, radius::Number)::Subset{Mode}
    currentspace::RealSpace = correlations.outspace |> getcrystal |> getspace |> orthogonalspace
    physicalnorm = m -> lineartransform(currentspace, m |> getpos) |> norm
    return filter(p -> physicalnorm(p - center) < radius, physicalmodes)
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

function localregioninspection(center::Offset, physicalmodes::Subset{Mode}, radius::Number)::Tuple{Subset{Offset},FockSpace}
    seedingmodes::Subset{Mode} = circularregionmodes(center, physicalmodes, radius)
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

triangular = RealSpace([sqrt(3)/2 -1/2; 0. 1.]')
kspace = convert(MomentumSpace, triangular)

pa = triangular & [1/3, 2/3]
pb = triangular & [2/3, 1/3]
pc = (pa + pb) / 2
spatialsnappingcalibration((pa, pb, pc))

c6 = pointgrouptransform([cos(π/3) -sin(π/3); sin(π/3) cos(π/3)])
c3 = c6^2

unitcell = Subset(pa, pb)
crystal = Crystal(unitcell, [48, 48])
reciprocalhashcalibration(crystal.sizes)

modes::Subset{Mode} = quantize(:pos, unitcell, 1)
m0, m1 = members(modes)

tₙ = ComplexF64(-1.)
bonds::FockMap = bondmap([
    (m0, m1) => tₙ,
    (m0, m1 |> setattr(:offset => triangular & [-1, 0])) => tₙ,
    (m0, m1 |> setattr(:offset => triangular & [0, 1])) => tₙ])

energyspectrum = computeenergyspectrum(bonds, crystal=crystal)
energyspectrum |> visualize

groundstates::CrystalSpectrum = groundstatespectrum(energyspectrum, perunitcellfillings=1)
groundstates |> visualize
groundstateprojector = groundstates |> crystalprojector

C = idmap(groundstateprojector.outspace) - groundstateprojector

correlations = C

crystalfock = correlations.outspace

scale = Scale([2 0; 0 2])
blockresult = blocking(:action => scale, :correlations => correlations, :crystal => crystalfock |> getcrystal)

blockedcrystal::Crystal = blockresult[:crystal]
blockedcorrelations::FockMap = blockresult[:correlations]

crystalpoints::Subset{Offset} = latticepoints(blockedcrystal)
blockedmodes::Subset{Mode} = quantize(:pos, blockedcrystal.unitcell, 1)
physicalmodes::Subset{Mode} = spanoffset(blockedmodes, crystalpoints)
scaledtriangular = scale*triangular 

localregion,localfock = localregioninspection(scaledtriangular&[0,0] , physicalmodes, 2)
visualize(localregion)

localcorrelation = regioncorrelations(blockedcorrelations,localfock)
localeigspec = localcorrelation |> eigspech 
visualize(localeigspec)

function locaclRG(center)::Tuple{FockMap,FockMap}
    localregion,localfock = localregioninspection(scaledtriangular&center , physicalmodes, 2)
    visualize(localregion)

    localcorrelation = regioncorrelations(blockedcorrelations,localfock)
    localeigspec = localcorrelation |> eigspech 
    visualize(localeigspec)

    localmodesdict = localmodesgrouping(localcorrelation, 0.01)

    c6recenter = c6 |> recenter(scaledtriangular&center)
    c3recenter = c3 |> recenter(scaledtriangular&center)

    # finding seeds for local frozen (6 modes at the center)
    frozenseedingcenter::Offset = (blockedmodes |> getspace) & center
    frozenseedingregion,frozenseedingfock = localregioninspection(frozenseedingcenter, physicalmodes, 1)
    visualize(frozenseedingregion, title="Frozen Seeding Region", visualspace=euclidean(RealSpace, 2))

    regioncorrelations(blockedcorrelations,frozenseedingfock) |> eigspech |>visualize

    frozenseedsorig = reduce(+ ,FockMap(v,inspace = v |> getinspace |> orderedmodes |> setattr(:dumind => n) |> FockSpace,performpermute = false) for (n, (_ , v)) in enumerate(_findlocalspstates(statecorrelations = blockedcorrelations, regionfock = frozenseedingfock, 
    spectrumextractpredicate = v -> true, symmetry = identitytransform(2))))

    # Summing all the frozen seed at center and corners to form a local frozen seeds
    # frozenseeds = frozenseedsorig+frozenseedAs+frozenseedBs

    frozenseeds = idmap(localfock, localfock)[:,frozenseedsorig |> getoutspace] * frozenseedsorig

    # Wannierization for frozen
    wannierizedfrozens = localwannierization(localmodesdict[:frozen], frozenseeds)

    # finding seeds for local courier (6 modes at the corners)
    courierseedingcenterA::Offset = (blockedmodes |> getspace) & (center+[2/3, 1/3])
    courierseedingregionA,courierseedingfockA = localregioninspection(courierseedingcenterA, physicalmodes, 0.8)
    visualize(courierseedingregionA, title="Frozen Seeding Region A", visualspace=euclidean(RealSpace, 2))
    
    courierseedingfockB = FockSpace{Region}(m for m in c6recenter*courierseedingfockA |> getoutspace)

    courierseedingfockArot1 = FockSpace{Region}(m for m in c3recenter*courierseedingfockA |> getoutspace)
    courierseedingfockBrot1 = FockSpace{Region}(m for m in c6recenter*c3recenter*courierseedingfockA |> getoutspace)

    courierseedingfockArot2 = FockSpace{Region}(m for m in c3recenter*c3recenter*courierseedingfockA |> getoutspace)
    courierseedingfockBrot2 = FockSpace{Region}(m for m in c6recenter*c3recenter*c3recenter*courierseedingfockA |> getoutspace)

    #regioncorrelations(blockedcorrelations,frozenseedingfockA) |> eigspech |>visualize
    courierseedA = reduce(+ ,FockMap(v,inspace = v |> getinspace |> orderedmodes |> setattr(:dumind => n) |> FockSpace,performpermute = false) for (n, (_ , v)) in enumerate(_findlocalspstates(statecorrelations = blockedcorrelations, regionfock = courierseedingfockA, 
        spectrumextractpredicate = v -> v < 0.6, symmetry = identitytransform(2))))

    courierseedB = reduce(+ ,FockMap(v,inspace = v |> getinspace |> orderedmodes |> setattr(:dumind => n) |> FockSpace,performpermute = false) for (n, (_ , v)) in enumerate(_findlocalspstates(statecorrelations = blockedcorrelations, regionfock = courierseedingfockB, 
        spectrumextractpredicate = v -> 0.4 < v, symmetry = identitytransform(2))))

    courierseedArot1 = reduce(+ ,FockMap(v,inspace = v |> getinspace |> orderedmodes |> setattr(:dumind => n) |> FockSpace,performpermute = false) for (n, (_ , v)) in enumerate(_findlocalspstates(statecorrelations = blockedcorrelations, regionfock = courierseedingfockArot1, 
        spectrumextractpredicate = v -> v < 0.6, symmetry = identitytransform(2))))

    courierseedBrot1 = reduce(+ ,FockMap(v,inspace = v |> getinspace |> orderedmodes |> setattr(:dumind => n) |> FockSpace,performpermute = false) for (n, (_ , v)) in enumerate(_findlocalspstates(statecorrelations = blockedcorrelations, regionfock = courierseedingfockBrot1, 
        spectrumextractpredicate = v -> 0.4 < v, symmetry = identitytransform(2))))

    courierseedArot2 = reduce(+ ,FockMap(v,inspace = v |> getinspace |> orderedmodes |> setattr(:dumind => n) |> FockSpace,performpermute = false) for (n, (_ , v)) in enumerate(_findlocalspstates(statecorrelations = blockedcorrelations, regionfock = courierseedingfockArot2, 
        spectrumextractpredicate = v -> v < 0.6, symmetry = identitytransform(2))))

    courierseedBrot2 = reduce(+ ,FockMap(v,inspace = v |> getinspace |> orderedmodes |> setattr(:dumind => n) |> FockSpace,performpermute = false) for (n, (_ , v)) in enumerate(_findlocalspstates(statecorrelations = blockedcorrelations, regionfock = courierseedingfockBrot2, 
        spectrumextractpredicate = v -> 0.4 < v, symmetry = identitytransform(2))))

    # Summing up all the symmetry related for corners
    courierseedAs = courierseedA + courierseedArot1 + courierseedArot2
    courierseedBs = courierseedB + courierseedBrot1 + courierseedBrot2

    courierseeds = courierseedAs+courierseedBs

    # Wannierization for courier
    wannierizedcouriers = _localwannierization(localmodesdict[:courier], courierseeds)

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

transformedphysicalmodes = transformedblockedcorrelationsRS |> getoutspace |> orderedmodes

trsasnformedlocalregion, trsasnformedlocalfock = localregioninspection(scaledtriangular&[0,0],transformedphysicalmodes, 2)
visualize(trsasnformedlocalregion, title="l", visualspace=euclidean(RealSpace, 2))

transformedlocaleigspec = transformedblockedcorrelationsRS[trsasnformedlocalfock, trsasnformedlocalfock] |> eigspech 

visualize(transformedlocaleigspec)
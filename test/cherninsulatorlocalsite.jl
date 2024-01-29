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
        return inmode |> setattr(:r => offset) |> setattr(:b => basis) |> setattr(:ind => ind)
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

modebydist = groupmodesbydist(region = localregion,regionfock = localfock,center = center)

function locaclRG(center_pt)::Tuple{FockMap,FockMap}
    center = Point(center_pt, scaledtriangular)

    localregion,localfock = localregioninspection(center , physicalmodes, 2,blockedcrystal)
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
    frozenseedingregion,frozenseedingfock = localregioninspection(frozenseedingcenter, physicalmodes, 1, blockedcrystal)
    visualize(frozenseedingregion, title="Frozen Seeding Region", visualspace=euclidean(RealSpace, 2))

    regioncorrelations(blockedcorrelations,frozenseedingfock) |> eigspech |>visualize

    frozenseedsorig = reduce(+ ,FockMap(v,inspace = v |> getinspace |> orderedmodes |> setattr(:dumind => n) |> FockSpace,performpermute = false) for (n, (_ , v)) in enumerate(_findlocalspstates(statecorrelations = blockedcorrelations, regionfock = frozenseedingfock, 
    spectrumextractpredicate = v -> true, symmetry = identitytransform(2))))

    extendedfrozenseeds = idmap(localfock, localfock)[:,frozenseedsorig |> getoutspace] * frozenseedsorig

    wannierizedfrozens = _localwannierization(localmodesdict[:frozen], extendedfrozenseeds)

    # finding seeds for local courier (6 modes at the corners)
    courierseedingcenterA::Offset = Point(center_pt+[2/3, 1/3] ,  (blockedmodes |> getspace))
    courierseedingregionA,courierseedingfockA = localregioninspection(courierseedingcenterA, physicalmodes, 0.1, blockedcrystal)
    visualize(courierseedingregionA, title="Courier Seeding Region A", visualspace=euclidean(RealSpace, 2))
    courierseedingfockB = FockSpace{Region}(m for m in c6recenter*courierseedingfockA |> getoutspace)


    regioncorrelations(blockedcorrelations,courierseedingfockB) |> eigspech |>visualize

    # for A 
    courierseedA =  findlocalspstates(statecorrelations = blockedcorrelations, regionfock = courierseedingfockA, 
        spectrumextractpredicate = v -> true, symmetry = identitytransform(2))[1]
    # for B 
    courierseedB =  findlocalspstates(statecorrelations = blockedcorrelations, regionfock = courierseedingfockB, 
        spectrumextractpredicate = v -> true, symmetry = identitytransform(2))[1]

    # Summing up all the symmetry related for corners
    courierseedAs = courierseedA + (c3recenter*courierseedA.outspace)*courierseedA*(c3recenter*courierseedA.inspace)' + (c3recenter^2*courierseedA.outspace)*courierseedA*(c3recenter^2*courierseedA.inspace)'
    courierseedBs = courierseedB + (c3recenter*courierseedB.outspace)*courierseedB*(c3recenter*courierseedB.inspace)' + (c3recenter^2*courierseedB.outspace)*courierseedB*(c3recenter^2*courierseedB.inspace)'

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
    courierseedAprimes = courierseedA1s + courierseedA2s

    # Summing up all the symmetry related for B (boundary)
    courierseedB1s = courierseedB1 + (c3recenter*courierseedB1.outspace)*courierseedB1*(c3recenter*courierseedB1.inspace)' + (c3recenter^2*courierseedB1.outspace)*courierseedB1*(c3recenter^2*courierseedB1.inspace)'
    courierseedB2s = courierseedB2 + (c3recenter*courierseedB2.outspace)*courierseedB2*(c3recenter*courierseedB2.inspace)' + (c3recenter^2*courierseedB2.outspace)*courierseedB2*(c3recenter^2*courierseedB2.inspace)'
    courierseedBprimes = courierseedB1s + courierseedB2s

    courierseeds = courierseedAs + courierseedBs + courierseedAprimes + courierseedBprimes

    extendedcourierseeds = idmap(localfock, localfock)[:, courierseeds |> getoutspace] * courierseeds

    wannierizedcouriers = _localwannierization(localmodesdict[:courier], extendedcourierseeds)

    return wannierizedfrozens, wannierizedcouriers
end

center_pt1 = [-2,-2]
wannierizedfrozens1, wannierizedcouriers1 = locaclRG(center_pt1)

center_pt2 = [-1,0]
wannierizedfrozens2, wannierizedcouriers2 = locaclRG(center_pt2)

center_pt3 = [0,2]
wannierizedfrozens3, wannierizedcouriers3 = locaclRG(center_pt3)

center_pt4 = [1,1]
wannierizedfrozens4, wannierizedcouriers4 = locaclRG(center_pt4)

center_pt5 = [2,0]
wannierizedfrozens5, wannierizedcouriers5 = locaclRG(center_pt5)

center_pt6 = [0,-1]
wannierizedfrozens6, wannierizedcouriers6 = locaclRG(center_pt6)

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

# transformedphysicalmodes = transformedblockedcorrelationsRS |> getoutspace |> orderedmodes


# trsasnformedlocalregion, trsasnformedlocalfock = localregioninspection(Point([0,0], scaledtriangular),transformedphysicalmodes, 2, blockedcrystal)
# visualize(trsasnformedlocalregion, title="l", visualspace=euclidean(RealSpace, 2))

# transformedlocaleigspec = transformedblockedcorrelationsRS[trsasnformedlocalfock, trsasnformedlocalfock] |> eigspech 

# visualize(transformedlocaleigspec)

function locaclRGsecond(center_pt)::Tuple{FockMap,FockMap}
    center = Point(center_pt, scaledtriangular)
    # inspecting local region spectrum
    transformedphysicalmodes = transformedblockedcorrelationsRS |> getoutspace |> orderedmodes
    
    trsasnformedlocalregion, trsasnformedlocalfock = localregioninspection(center ,transformedphysicalmodes, 2, blockedcrystal)
    visualize(trsasnformedlocalregion, title="l", visualspace=euclidean(RealSpace, 2))
    
    transformedlocaleigspec = transformedblockedcorrelationsRS[trsasnformedlocalfock, trsasnformedlocalfock] |> eigspech 
    visualize(transformedlocaleigspec)
    
    
    transformedlocalmodesdict = localmodesgrouping(transformedblockedcorrelationsRS[trsasnformedlocalfock, trsasnformedlocalfock], 0.01)
    
    # c6recenter = c6 |> recenter(scaledtriangular&center)
    # c3recenter = c3 |> recenter(scaledtriangular&center)
    
    # finding seeds for local frozen (6 modes at the centers)
    frozenseedingcenter::Offset = Point(center_pt, scaledtriangular)
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
    
    
    # finding seeds for local courier (6 modes closer to the center)
    courierseedingcenternA1::Offset = Point(center_pt+[2/3, 1/3] ,  scaledtriangular)
    courierseedingregionnA1RG,courierseedingfocknA1RG = localregioninspection(courierseedingcenternA1, transformedphysicalmodes, 0.1, blockedcrystal)
    visualize(courierseedingregionnA1RG, title="Courier Seeding Region A", visualspace=euclidean(RealSpace, 2))
    
    courierseedingcenternB1::Offset = Point(center_pt+[1/3, 2/3] ,  scaledtriangular)
    courierseedingregionnB1RG,courierseedingfocknB1RG = localregioninspection(courierseedingcenternB1, transformedphysicalmodes, 0.1, blockedcrystal)
    visualize(courierseedingregionnB1RG, title="Courier Seeding Region B", visualspace=euclidean(RealSpace, 2))
    
    # for A 
    courierseednA1 = reduce(+ ,FockMap(v,inspace = v |> getinspace |> orderedmodes |> setattr(:dumind => n) |> FockSpace,performpermute = false) for (n, (_ , v)) in enumerate(findlocalseeds(statecorrelations = transformedblockedcorrelationsRS, regionfock = courierseedingfocknA1RG, 
    spectrumextractpredicate = v -> true, symmetry = identitytransform(2))))
    
    # for B 
    courierseednB1 =  reduce(+ ,FockMap(v,inspace = v |> getinspace |> orderedmodes |> setattr(:dumind => n) |> FockSpace,performpermute = false) for (n, (_ , v)) in enumerate(findlocalseeds(statecorrelations = transformedblockedcorrelationsRS, regionfock = courierseedingfocknB1RG, 
    spectrumextractpredicate = v -> true, symmetry = identitytransform(2))))
    
    
    courierseedingcenternA2::Offset = Point(center_pt+[-1/3, -2/3] ,  scaledtriangular)
    courierseedingregionnA2RG,courierseedingfocknA2RG = localregioninspection(courierseedingcenternA2, transformedphysicalmodes, 0.1, blockedcrystal)
    visualize(courierseedingregionnA2RG, title="Courier Seeding Region A", visualspace=euclidean(RealSpace, 2))
    
    courierseedingcenternB2::Offset = Point(center_pt+[1/3, -1/3] ,  scaledtriangular)
    courierseedingregionnB2RG,courierseedingfocknB2RG = localregioninspection(courierseedingcenternB2, transformedphysicalmodes, 0.1, blockedcrystal)
    visualize(courierseedingregionnB2RG, title="Courier Seeding Region B", visualspace=euclidean(RealSpace, 2))
    
    # for A 
    courierseednA2 =  reduce(+ ,FockMap(v,inspace = v |> getinspace |> orderedmodes |> setattr(:dumind => n) |> FockSpace,performpermute = false) for (n, (_ , v)) in enumerate(findlocalseeds(statecorrelations = transformedblockedcorrelationsRS, regionfock = courierseedingfocknA2RG, 
    spectrumextractpredicate = v -> true, symmetry = identitytransform(2))))
    # for B 
    courierseednB2 =  reduce(+ ,FockMap(v,inspace = v |> getinspace |> orderedmodes |> setattr(:dumind => n) |> FockSpace,performpermute = false) for (n, (_ , v)) in enumerate(findlocalseeds(statecorrelations = transformedblockedcorrelationsRS, regionfock = courierseedingfocknB2RG, 
    spectrumextractpredicate = v -> true, symmetry = identitytransform(2))))
    
    
    courierseedingcenternA3::Offset = Point(center_pt+[-1/3, 1/3] ,  scaledtriangular)
    courierseedingregionnA3RG,courierseedingfocknA3RG = localregioninspection(courierseedingcenternA3, transformedphysicalmodes, 0.1, blockedcrystal)
    visualize(courierseedingregionnA3RG, title="Courier Seeding Region A", visualspace=euclidean(RealSpace, 2))
    
    courierseedingcenternB3::Offset = Point(center_pt+[-2/3, -1/3] ,  scaledtriangular)
    courierseedingregionnB3RG,courierseedingfocknB3RG = localregioninspection(courierseedingcenternB3, transformedphysicalmodes, 0.1, blockedcrystal)
    visualize(courierseedingregionnB3RG, title="Courier Seeding Region B", visualspace=euclidean(RealSpace, 2))
    
    # for A 
    courierseednA3 =  reduce(+ ,FockMap(v,inspace = v |> getinspace |> orderedmodes |> setattr(:dumind => n) |> FockSpace,performpermute = false) for (n, (_ , v)) in enumerate(findlocalseeds(statecorrelations = transformedblockedcorrelationsRS, regionfock = courierseedingfocknA3RG, 
    spectrumextractpredicate = v -> true, symmetry = identitytransform(2))))
    # for B 
    courierseednB3 =  reduce(+ ,FockMap(v,inspace = v |> getinspace |> orderedmodes |> setattr(:dumind => n) |> FockSpace,performpermute = false) for (n, (_ , v)) in enumerate(findlocalseeds(statecorrelations = transformedblockedcorrelationsRS, regionfock = courierseedingfocknB3RG, 
    spectrumextractpredicate = v -> true, symmetry = identitytransform(2))))
    
    # finding seeds for local courier (6 modes at the corners)
    courierseedingcenterA1::Offset = Point(center_pt+[0 + 1/6, -3/4 + 1/12], scaledtriangular)
    courierregionA1RG, courierseedingfockA1RG = localregioninspection(courierseedingcenterA1, transformedphysicalmodes, 0.1, blockedcrystal)
    visualize(courierregionA1RG, visualspace=euclidean(RealSpace, 2))
    
    courierseedingcenterB1::Offset = Point(center_pt+[0 - 1/6, -3/4 - 1/12], scaledtriangular)
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
    
    courierseedingcenterB2::Offset = Point(center_pt+[-3/4 + 1/12, 0 + 1/6], scaledtriangular)
    courierregionB2RG, courierseedingfockB2RG = localregioninspection(courierseedingcenterB2, transformedphysicalmodes, 0.1, blockedcrystal)
    courierB2RGspec = transformedblockedcorrelationsRS[courierseedingfockB2RG, courierseedingfockB2RG] |> eigspech 
    visualize(courierB2RGspec)
    
    courierseedingcenterA2::Offset = Point(center_pt+[-3/4 - 1/12, 0 - 1/6], scaledtriangular)
    courierregionA2RG, courierseedingfockA2RG = localregioninspection(courierseedingcenterA2, transformedphysicalmodes, 0.1, blockedcrystal)
    courierA2RGspec = transformedblockedcorrelationsRS[courierseedingfockA2RG, courierseedingfockA2RG] |> eigspech 
    visualize(courierA2RGspec)
    
    courierseedA2RG = reduce(+ ,FockMap(v,inspace = v |> getinspace |> orderedmodes |> setattr(:dumind => n) |> FockSpace,performpermute = false) for (n, (_ , v)) in enumerate(findlocalseeds(statecorrelations = transformedblockedcorrelationsRS, regionfock = courierseedingfockA2RG, 
    spectrumextractpredicate = v -> true, symmetry = identitytransform(2))))
    
    courierseedB2RG = reduce(+ ,FockMap(v,inspace = v |> getinspace |> orderedmodes |> setattr(:dumind => n) |> FockSpace,performpermute = false) for (n, (_ , v)) in enumerate(findlocalseeds(statecorrelations = transformedblockedcorrelationsRS, regionfock = courierseedingfockB2RG, 
    spectrumextractpredicate = v -> true, symmetry = identitytransform(2))))
    
    courierseedingcenterB3::Offset = Point(center_pt+[3/4 + 1/12, 3/4 - 1/12], scaledtriangular)
    courierregionB3RG, courierseedingfockB3RG = localregioninspection(courierseedingcenterB3, transformedphysicalmodes, 0.1, blockedcrystal)
    courierB3RGspec = transformedblockedcorrelationsRS[courierseedingfockB3RG, courierseedingfockB3RG] |> eigspech 
    visualize(courierB3RGspec)
    
    courierseedingcenterA3::Offset = Point(center_pt+[3/4 - 1/12, 3/4 + 1/12], scaledtriangular)
    courierregionA3RG, courierseedingfockA3RG = localregioninspection(courierseedingcenterA3, transformedphysicalmodes, 0.1, blockedcrystal)
    courierA3RGspec = transformedblockedcorrelationsRS[courierseedingfockA3RG, courierseedingfockA3RG] |> eigspech 
    visualize(courierA3RGspec)
    
    courierseedA3RG = reduce(+ ,FockMap(v,inspace = v |> getinspace |> orderedmodes |> setattr(:dumind => n) |> FockSpace,performpermute = false) for (n, (_ , v)) in enumerate(findlocalseeds(statecorrelations = transformedblockedcorrelationsRS, regionfock = courierseedingfockA3RG, 
    spectrumextractpredicate = v -> true, symmetry = identitytransform(2))))
    
    courierseedB3RG = reduce(+ ,FockMap(v,inspace = v |> getinspace |> orderedmodes |> setattr(:dumind => n) |> FockSpace,performpermute = false) for (n, (_ , v)) in enumerate(findlocalseeds(statecorrelations = transformedblockedcorrelationsRS, regionfock = courierseedingfockB3RG, 
    spectrumextractpredicate = v -> true, symmetry = identitytransform(2))))
    
    
    courierseedRG = courierseednA1 + courierseednB1 + courierseednA2 + courierseednB2 + courierseednA3+ courierseednB3 + courierseedA1RG + courierseedB1RG + courierseedA2RG + courierseedB2RG + courierseedA3RG + courierseedB3RG
    
    extendedcourierseedsRG = idmap(trsasnformedlocalfock, trsasnformedlocalfock)[:,courierseedRG |> getoutspace] * courierseedRG
    
    # Wannierization for couriers
    wannierizedcouriersRG = _localwannierization(transformedlocalmodesdict[:courier], extendedcourierseedsRG)
    
    return wannierizedfrozensRG, wannierizedcouriersRG
end

transformedblockedcorrelationsRS |> getinspace 

center_pt = [1,0]
center = Point(center_pt, scaledtriangular)
# inspecting local region spectrum
transformedphysicalmodes = transformedblockedcorrelationsRS |> getoutspace |> orderedmodes
    
trsasnformedlocalregion, trsasnformedlocalfock = localregioninspection(center ,transformedphysicalmodes, 2, blockedcrystal)
visualize(trsasnformedlocalregion, title="l", visualspace=euclidean(RealSpace, 2))
    
transformedlocaleigspec = transformedblockedcorrelationsRS[trsasnformedlocalfock, trsasnformedlocalfock] |> eigspech 
visualize(transformedlocaleigspec)

center_ptA = [1,0]
wannierizedfrozensRGA, wannierizedcouriersRGA = locaclRGsecond(center_ptA)::Tuple{FockMap,FockMap}

wannierizedcouriersRGA[:,3] |> Zipper.columnspec |> visualize

center_ptB = [0,1]
wannierizedfrozensRGB, wannierizedcouriersRGB = locaclRGsecond(center_ptB)::Tuple{FockMap,FockMap}

center_ptC = [-1,-1]
wannierizedfrozensRGC, wannierizedcouriersRGC = locaclRGsecond(center_ptC)::Tuple{FockMap,FockMap}


localunitaryRGsecondA = wannierizedfrozensRGA + wannierizedcouriersRGA
localunitaryRGsecondB = wannierizedfrozensRGB + wannierizedcouriersRGB
localunitaryRGsecondC = wannierizedfrozensRGC + wannierizedcouriersRGC

extendedisosecond = ((((transformedblockedcorrelationsRS|> getoutspace) - (localunitaryRGsecondA  |> getoutspace) - (localunitaryRGsecondB  |> getoutspace) - (localunitaryRGsecondC  |> getoutspace)) |> idmap)
                 + wannierizedcouriersRGA + wannierizedcouriersRGB + wannierizedcouriersRGC)

transformedblockedcorrelationsRSsecond = extendedisosecond'*transformedblockedcorrelationsRS*extendedisosecond

# inspecting local region spectrum
center_pt = [0,0]
transformedphysicalmodessecond = transformedblockedcorrelationsRSsecond |> getoutspace |> orderedmodes

trsasnformedlocalregionsecond, trsasnformedlocalfocksecond = localregioninspection(Point([0,0], scaledtriangular),transformedphysicalmodessecond, 2, blockedcrystal)
visualize(trsasnformedlocalregionsecond, title="l", visualspace=euclidean(RealSpace, 2))

transformedlocaleigspec = transformedblockedcorrelationsRSsecond[trsasnformedlocalfocksecond, trsasnformedlocalfocksecond] |> eigspech 
visualize(transformedlocaleigspec)

transformedlocalmodesdictsecond = localmodesgrouping(transformedblockedcorrelationsRSsecond[trsasnformedlocalfocksecond, trsasnformedlocalfocksecond], 0.003)


# visualize((wannierizedfrozensRGA'*transformedblockedcorrelationsRS[trsasnformedlocalfock, trsasnformedlocalfock]*wannierizedfrozensRGA) |> eigspech )

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
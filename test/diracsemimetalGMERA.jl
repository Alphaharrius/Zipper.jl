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

function translatelocalfockmap(;
    localfockmap::FockMap,
    offset::Offset)::FockMap

    outspace = FockSpace((mode+offset) for mode in localfockmap |> getoutspace)
    inspace = FockSpace((mode+offset) for mode in localfockmap |> getinspace)
    representation = (localfockmap |> rep)

    return FockMap(outspace, inspace, representation)
end

function translatedlocalfockmaplist(;
    center::Offset,
    localfockmap::FockMap,
    transceterlist::Set{Offset})::Vector{FockMap}
    localfockmaplist::Vector{FockMap} = []
    for transcenter in transceterlist
        shift::Offset = transcenter - center
        # println(transcenter)
        # println(shift)
        push!(localfockmaplist,translatelocalfockmap(localfockmap = localfockmap,offset = shift))
    end
    return localfockmaplist
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
crystal = Crystal(unitcell, [96, 96])
reciprocalhashcalibration(crystal.sizes)

modes::Subset{Mode} = quantize(:pos, unitcell, 1)
m0, m1 = members(modes)

tₙ = ComplexF64(-1.)
bonds::FockMap = bondmap([
    (m0, m1) => tₙ,
    (m0, setattr(m1, :offset => Point([-1, 0], triangular))) => tₙ,
    (m0, setattr(m1, :offset => Point([0, 1], triangular))) => tₙ])

energyspectrum = computeenergyspectrum(bonds, crystal=crystal)
energyspectrum |> visualize

groundstates::CrystalSpectrum = groundstatespectrum(energyspectrum, perunitcellfillings=1)
groundstates |> visualize
groundstateprojector = groundstates |> crystalprojector

C = idmap(groundstateprojector.outspace) - groundstateprojector

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

firstlayer = Set([Point([-2,-2],scaledtriangular), Point([-1,0],scaledtriangular), Point([0,2],scaledtriangular), Point([0,-1],scaledtriangular), Point([1,1],scaledtriangular), Point([2,0],scaledtriangular)])
secondlayer = Set([Point([-1,-1],scaledtriangular), Point([0,1],scaledtriangular), Point([1,0],scaledtriangular)])
thirdlayer = Set([Point([0,0],scaledtriangular)])

# firstlayer = Set([[-2,-2], [-1,0], [0,2], [0,-1], [1,1], [2,0]])
# secondlayer = Set([[-1,-1], [0,1], [1,0]])
# thirdlayer = Set([[0,0]])

firstlayers = Set([Point([-2,-2],scaledtriangular), Point([-1,0],scaledtriangular), Point([0,2],scaledtriangular), Point([0,-1],scaledtriangular), Point([1,1],scaledtriangular), Point([2,0],scaledtriangular)])
secondlayers = Set([Point([-1,-1],scaledtriangular), Point([0,1],scaledtriangular), Point([1,0],scaledtriangular)])
thirdlayers = Set([Point([0,0],scaledtriangular)])

# firstlayers = Set([[-2,-2], [-1,0], [0,2], [0,-1], [1,1], [2,0]])
# secondlayers = Set([[-1,-1], [0,1], [1,0]])
# thirdlayers = Set([[0,0]])

for i in range(-6,6)
    for j in range(-6,6)
        firstlayers = union(firstlayers,Set([r + Point(i*[1,2]+j*[2,1], scaledtriangular) for r in firstlayer]))
        secondlayers = union(secondlayers,Set([r + Point(i*[1,2]+j*[2,1], scaledtriangular) for r in secondlayer]))
        thirdlayers = union(thirdlayers,Set([r + Point(i*[1,2]+j*[2,1], scaledtriangular) for r in thirdlayer]))
    end
end

# for i in range(-7,7)
#     for j in range(-7,7)
#         firstlayers = union(firstlayers,Set([r+i*[1,2]+j*[2,1] for r in firstlayer]))
#         secondlayers = union(secondlayers,Set([r+i*[1,2]+j*[2,1] for r in secondlayer]))
#         thirdlayers = union(thirdlayers,Set([r+i*[1,2]+j*[2,1] for r in thirdlayer]))
#     end
# end

function locaclRGfirst(center::Offset,statecorrelation::FockMap, physicalmodes::Subset{Mode}, modegroupingthreshold::Float64, radius::Number)::Tuple{FockMap,FockMap}
    localregion,localfock = localregioninspection(center , physicalmodes, radius,blockedcrystal)
    visualize(localregion)

    localcorrelation = regioncorrelations(statecorrelation,localfock)
    localeigspec = localcorrelation |> eigspech 
    visualize(localeigspec)

    localmodesdict = localmodesgrouping(localcorrelation, modegroupingthreshold)
    modebydist = groupmodesbydist(region = localregion,regionfock = localfock,center = center)  

    # finding seeds for local frozen (6 modes at the center)
    frozenseedingfocks = (frozenseedingfock for (_,frozenseedingfock) in (localregioninspection(modewithdist[2] |> getpos, physicalmodes, 0.1, blockedcrystal) for modewithdist in modebydist[3]))
    frozenseeds = sum([reduce(+ ,FockMap(v,inspace = v |> getinspace |> orderedmodes |> setattr(:dumind => n) |> FockSpace,performpermute = false) for (n, (_ , v)) in enumerate(_findlocalspstates(statecorrelations = statecorrelation, regionfock = frozenseedingfock, 
    spectrumextractpredicate = v -> true, symmetry = identitytransform(2)))) for frozenseedingfock in frozenseedingfocks])

    extendedfrozenseeds = idmap(localfock, localfock)[:,frozenseeds |> getoutspace] * frozenseeds

    wannierizedfrozens = _localwannierization(localmodesdict[:frozen], extendedfrozenseeds)

    # finding seeds for local courier (18 (12+6) modes at "boundary")
    courierseedingfocks1 = (courierseedingfock for (_,courierseedingfock) in (localregioninspection(modewithdist[2] |> getpos, physicalmodes, 0.1, blockedcrystal) for modewithdist in modebydist[1]))
    courierseedingfocks2 = (courierseedingfock for (_,courierseedingfock) in (localregioninspection(modewithdist[2] |> getpos, physicalmodes, 0.1, blockedcrystal) for modewithdist in modebydist[2]))

    courierseeds1 = sum([reduce(+ ,FockMap(v,inspace = v |> getinspace |> orderedmodes |> setattr(:dumind => n) |> FockSpace,performpermute = false) for (n, (_ , v)) in enumerate(_findlocalspstates(statecorrelations = statecorrelation, regionfock = courierseedingfock, 
    spectrumextractpredicate = v -> true, symmetry = identitytransform(2)))) for courierseedingfock in courierseedingfocks1])
    courierseeds2 = sum([reduce(+ ,FockMap(v,inspace = v |> getinspace |> orderedmodes |> setattr(:dumind => n) |> FockSpace,performpermute = false) for (n, (_ , v)) in enumerate(_findlocalspstates(statecorrelations = statecorrelation, regionfock = courierseedingfock, 
    spectrumextractpredicate = v -> true, symmetry = identitytransform(2)))) for courierseedingfock in courierseedingfocks2])

    courierseeds = courierseeds1 + courierseeds2
    extendedcourierseeds = idmap(localfock, localfock)[:, courierseeds |> getoutspace] * courierseeds

    wannierizedcouriers = _localwannierization(localmodesdict[:courier], extendedcourierseeds)

    return wannierizedfrozens, wannierizedcouriers
end

function locaclRGsecond(center::Offset,statecorrelation::FockMap, physicalmodes::Subset{Mode}, modegroupingthreshold::Float64, radius::Number)::Tuple{FockMap,FockMap}
    localregion, localfock = localregioninspection(center ,physicalmodes, radius, blockedcrystal)
    visualize(localregion, title="l", visualspace=euclidean(RealSpace, 2))
    
    localeigspec = statecorrelation[localfock, localfock] |> eigspech 
    visualize(localeigspec)
    
    localmodesdict = localmodesgrouping(statecorrelation[localfock, localfock], modegroupingthreshold)
    modebydist = groupmodesbydist(region = localregion,regionfock = localfock,center = center)
    
    # finding seeds for local frozen (6 modes at the center)
    frozenseedingfocks = (frozenseedingfock for (_,frozenseedingfock) in (localregioninspection(modewithdist[2] |> getpos, physicalmodes, 0.1, blockedcrystal) for modewithdist in modebydist[3]))
    frozenseeds = sum([reduce(+ ,FockMap(v,inspace = v |> getinspace |> orderedmodes |> setattr(:dumind => n) |> FockSpace,performpermute = false) for (n, (_ , v)) in enumerate(findlocalseeds(statecorrelations = statecorrelation, regionfock = frozenseedingfock, 
    spectrumextractpredicate = v -> true, symmetry = identitytransform(2)))) for frozenseedingfock in frozenseedingfocks])

    extendedfrozenseeds = idmap(localfock, localfock)[:,frozenseeds |> getoutspace] * frozenseeds

    wannierizedfrozens = _localwannierization(localmodesdict[:frozen], extendedfrozenseeds)

    # finding seeds for local courier (18 (12+6) modes at "boundary")
    courierseedingfocks1 = (courierseedingfock for (_,courierseedingfock) in (localregioninspection(modewithdist[2] |> getpos, physicalmodes, 0.1, blockedcrystal) for modewithdist in modebydist[1]))
    courierseedingfocks2 = (courierseedingfock for (_,courierseedingfock) in (localregioninspection(modewithdist[2] |> getpos, physicalmodes, 0.1, blockedcrystal) for modewithdist in modebydist[2]))

    courierseeds1 = sum([reduce(+ ,FockMap(v,inspace = v |> getinspace |> orderedmodes |> setattr(:dumind => n) |> FockSpace,performpermute = false) for (n, (_ , v)) in enumerate(findlocalseeds(statecorrelations = statecorrelation, regionfock = courierseedingfock, 
    spectrumextractpredicate = v -> true, symmetry = identitytransform(2)))) for courierseedingfock in courierseedingfocks1])
    courierseeds2 = sum([reduce(+ ,FockMap(v,inspace = v |> getinspace |> orderedmodes |> setattr(:dumind => n) |> FockSpace,performpermute = false) for (n, (_ , v)) in enumerate(findlocalseeds(statecorrelations = statecorrelation, regionfock = courierseedingfock, 
    spectrumextractpredicate = v -> true, symmetry = identitytransform(2)))) for courierseedingfock in courierseedingfocks2])

    courierseeds = courierseeds1 + courierseeds2
    extendedcourierseeds = idmap(localfock, localfock)[:, courierseeds |> getoutspace] * courierseeds
    

    # Wannierization for couriers
    wannierizedcouriers = _localwannierization(localmodesdict[:courier], extendedcourierseeds)
    
    return wannierizedfrozens, wannierizedcouriers
end

function locaclRGthird(center::Offset,statecorrelation::FockMap, physicalmodes::Subset{Mode}, modegroupingthreshold::Float64, radius::Number)::Tuple{FockMap,FockMap}
    localregion, localfock = localregioninspection(center ,physicalmodes, radius, blockedcrystal)
    visualize(localregion, title="l", visualspace=euclidean(RealSpace, 2))
    
    localeigspec = statecorrelation[localfock, localfock] |> eigspech 
    visualize(localeigspec)
    
    localmodesdict = localmodesgrouping(statecorrelation[localfock, localfock], modegroupingthreshold)
    modebydist = groupmodesbydist(region = localregion,regionfock = localfock,center = center)
    
    # finding seeds for local frozen (6 modes at the center)
    frozenseedingfocks = (frozenseedingfock for (_,frozenseedingfock) in (localregioninspection(modewithdist[2] |> getpos, physicalmodes, 0.1, blockedcrystal) for modewithdist in modebydist[2]))
    frozenseeds = sum([reduce(+ ,FockMap(v,inspace = v |> getinspace |> orderedmodes |> setattr(:dumind => n) |> FockSpace,performpermute = false) for (n, (_ , v)) in enumerate(findlocalseeds(statecorrelations = statecorrelation, regionfock = frozenseedingfock, 
    spectrumextractpredicate = v -> true, symmetry = identitytransform(2)))) for frozenseedingfock in frozenseedingfocks])

    extendedfrozenseeds = idmap(localfock, localfock)[:,frozenseeds |> getoutspace] * frozenseeds

    wannierizedfrozens = _localwannierization(localmodesdict[:frozen], extendedfrozenseeds)

    # finding seeds for local courier (18 (12+6) modes at "boundary")
    courierseedingfocks = (courierseedingfock for (_,courierseedingfock) in (localregioninspection(modewithdist[2] |> getpos, physicalmodes, 0.1, blockedcrystal) for modewithdist in modebydist[1]))
    
    courierseeds = sum([reduce(+ ,FockMap(v,inspace = v |> getinspace |> orderedmodes |> setattr(:dumind => n) |> FockSpace,performpermute = false) for (n, (_ , v)) in enumerate(findlocalseeds(statecorrelations = statecorrelation, regionfock = courierseedingfock, 
    spectrumextractpredicate = v -> true, symmetry = identitytransform(2)))) for courierseedingfock in courierseedingfocks])
    
    extendedcourierseeds = idmap(localfock, localfock)[:, courierseeds |> getoutspace] * courierseeds
    

    # Wannierization for couriers
    wannierizedcouriers = _localwannierization(localmodesdict[:courier], extendedcourierseeds)
    
    return wannierizedfrozens, wannierizedcouriers
end

# inspecting local region spectrum
center_pt = [-2,2]
center = Point(center_pt, scaledtriangular)

localregion,localfock = localregioninspection(center , physicalmodes, 2,blockedcrystal)
visualize(localregion)

localcorrelation = regioncorrelations(blockedcorrelations,localfock)
localeigspec = localcorrelation |> eigspech 
visualize(localeigspec, title="Spectrum of Dirac semi-metal to start with")

#
center = Point([-2,-2],scaledtriangular)
refwannierfn = locaclRGfirst(center,blockedcorrelations, physicalmodes, 0.001, 2)
localunitaryref = refwannierfn[1]+refwannierfn[2]

localunitarylist = translatedlocalfockmaplist(center = center,localfockmap = localunitaryref,transceterlist=firstlayers)
wanniercourierlist = translatedlocalfockmaplist(center = center,localfockmap = refwannierfn[2],transceterlist=firstlayers)
localunitary = sum(localunitarylist)
wanniercouriers = sum(wanniercourierlist)

couriercorrelationsecond = wanniercouriers'*regioncorrelations(blockedcorrelations,localunitary |> getoutspace)*wanniercouriers

# inspecting local region spectrum
center = Point([-1,-1], scaledtriangular)
transformedphysicalmodesnext = couriercorrelationsecond |> getoutspace |> orderedmodes

trsasnformedlocalregionnext, trsasnformedlocalfocknext = localregioninspection(center,transformedphysicalmodesnext, 2, blockedcrystal)
visualize(trsasnformedlocalregionnext, title="l", visualspace=euclidean(RealSpace, 2))

localeigspecnext = couriercorrelationsecond[trsasnformedlocalfocknext, trsasnformedlocalfocknext] |> eigspech 
visualize(localeigspecnext)

#
physicalmodessecond = couriercorrelationsecond |> getoutspace |> orderedmodes

refwannierfnsecond = locaclRGsecond(center,couriercorrelationsecond, physicalmodessecond, 0.003, 2)
localunitaryrefsecond = refwannierfnsecond[1]+refwannierfnsecond[2]

localunitarylistsecond = translatedlocalfockmaplist(center = center,localfockmap = localunitaryrefsecond,transceterlist=secondlayers)
wanniercourierlistsecond = translatedlocalfockmaplist(center = center,localfockmap = refwannierfnsecond[2],transceterlist=secondlayers)
localunitarysecond = sum(localunitarylistsecond)
wanniercourierssecond = sum(wanniercourierlistsecond)

couriercorrelationthird = wanniercourierssecond'*couriercorrelationsecond[localunitarysecond |> getoutspace, localunitarysecond |> getoutspace]*wanniercourierssecond

# inspecting local region spectrum
center = Point([0,0], scaledtriangular)
transformedphysicalmodesnext = couriercorrelationthird |> getoutspace |> orderedmodes

trsasnformedlocalregionnext, trsasnformedlocalfocknext = localregioninspection(center,transformedphysicalmodesnext, 2, blockedcrystal)
visualize(trsasnformedlocalregionnext, title="l", visualspace=euclidean(RealSpace, 2))

localeigspecnext = couriercorrelationthird[trsasnformedlocalfocknext, trsasnformedlocalfocknext] |> eigspech 
visualize(localeigspecnext)

#
physicalmodesthird = couriercorrelationthird |> getoutspace |> orderedmodes

refwannierfnthird = locaclRGthird(center,couriercorrelationthird, physicalmodesthird, 0.002, 2)
localunitaryrefthird = refwannierfnthird[1]+refwannierfnthird[2]

localunitarylistthird = translatedlocalfockmaplist(center = center,localfockmap = localunitaryrefthird,transceterlist=thirdlayers)
wanniercourierlistthird = translatedlocalfockmaplist(center = center,localfockmap = refwannierfnthird[2],transceterlist=thirdlayers)
localunitarythird = sum(localunitarylistthird)
wanniercouriersthird = sum(wanniercourierlistthird)

couriercorrelationRG = wanniercouriersthird'*couriercorrelationthird[localunitarythird |> getoutspace, localunitarythird |> getoutspace]*wanniercouriersthird

# inspecting local region spectrum
center = Point([-4,-4], scaledtriangular)
physicalmodesRG = couriercorrelationRG |> getoutspace |> orderedmodes

trsasnformedlocalregionnext, trsasnformedlocalfocknext = localregioninspection(center,physicalmodesRG, 3.5, blockedcrystal)
visualize(trsasnformedlocalregionnext, title="l", visualspace=euclidean(RealSpace, 2))

localeigspecnext = couriercorrelationRG[trsasnformedlocalfocknext, trsasnformedlocalfocknext] |> eigspech 
visualize(localeigspecnext,title="Spectrum after 1 RG step for Dirac semi-metal")

firstlayerRG = Set([Point([-4,-4],scaledtriangular), Point([-2,0],scaledtriangular), Point([0,4],scaledtriangular), Point([0,-2],scaledtriangular), Point([2,2],scaledtriangular), Point([4,0],scaledtriangular)])
secondlayerRG = Set([Point([-2,-2],scaledtriangular), Point([0,2],scaledtriangular), Point([2,0],scaledtriangular)])
thirdlayerRG = Set([Point([0,0],scaledtriangular)])

firstlayerRGs = Set([Point([-4,-4],scaledtriangular), Point([-2,0],scaledtriangular), Point([0,4],scaledtriangular), Point([0,-2],scaledtriangular), Point([2,2],scaledtriangular), Point([4,0],scaledtriangular)])
secondlayerRGs = Set([Point([-2,-2],scaledtriangular), Point([0,2],scaledtriangular), Point([2,0],scaledtriangular)])
thirdlayerRGs = Set([Point([0,0],scaledtriangular)])


for i in range(-1,1)
    for j in range(-1,1)
        firstlayerRGs = union(firstlayerRGs,Set([r + Point(i*[2,4]+j*[4,2], scaledtriangular) for r in firstlayerRG]))
        secondlayerRGs = union(secondlayerRGs,Set([r + Point(i*[2,4]+j*[4,2], scaledtriangular) for r in secondlayerRG]))
        thirdlayerRGs = union(thirdlayerRGs,Set([r + Point(i*[2,4]+j*[4,2], scaledtriangular) for r in thirdlayerRG]))
    end
end

#
physicalmodesRG = couriercorrelationRG |> getoutspace |> orderedmodes

center = Point([-4,-4],scaledtriangular)
refwannierfnRG = locaclRGsecond(center,couriercorrelationRG, physicalmodesRG , 0.003, 4)
localunitaryrefRG = refwannierfnRG[1]+refwannierfnRG[2]

localunitarylistRG = translatedlocalfockmaplist(center = center,localfockmap = localunitaryrefRG,transceterlist=firstlayerRGs)
wanniercourierlistRG = translatedlocalfockmaplist(center = center,localfockmap = refwannierfnRG[2],transceterlist=firstlayerRGs)
localunitaryRG = sum(localunitarylistRG)
wanniercouriersRG = sum(wanniercourierlistRG)

couriercorrelationRGsecond = wanniercouriersRG'*couriercorrelationRG[localunitaryRG |> getoutspace, localunitaryRG |> getoutspace]*wanniercouriersRG


# inspecting local region spectrum
center = Point([2,0], scaledtriangular)
transformedphysicalmodesnext = couriercorrelationRGsecond |> getoutspace |> orderedmodes

trsasnformedlocalregionnext, trsasnformedlocalfocknext = localregioninspection(center,transformedphysicalmodesnext, 4, blockedcrystal)
visualize(trsasnformedlocalregionnext, title="l", visualspace=euclidean(RealSpace, 2))

localeigspecnext = couriercorrelationRGsecond[trsasnformedlocalfocknext, trsasnformedlocalfocknext] |> eigspech 
visualize(localeigspecnext)

#
physicalmodesRGsecond = couriercorrelationRGsecond |> getoutspace |> orderedmodes

refwannierfnsecondRG = locaclRGsecond(center,couriercorrelationRGsecond, physicalmodesRGsecond, 0.005, 4)
localunitaryrefsecondRG = refwannierfnsecondRG[1]+refwannierfnsecondRG[2]

localunitarylistsecondRG = translatedlocalfockmaplist(center = center,localfockmap = localunitaryrefsecondRG,transceterlist=secondlayerRGs)
wanniercourierlistsecondRG = translatedlocalfockmaplist(center = center,localfockmap = refwannierfnsecondRG[2],transceterlist=secondlayerRGs)
localunitarysecondRG = sum(localunitarylistsecondRG)
wanniercourierssecondRG = sum(wanniercourierlistsecondRG)

couriercorrelationRGthird = wanniercourierssecondRG'*couriercorrelationRGsecond[localunitarysecondRG |> getoutspace, localunitarysecondRG |> getoutspace]*wanniercourierssecondRG

# inspecting local region spectrum
center = Point([0,0], scaledtriangular)
transformedphysicalmodesnext = couriercorrelationRGthird |> getoutspace |> orderedmodes

trsasnformedlocalregionnext, trsasnformedlocalfocknext = localregioninspection(center,transformedphysicalmodesnext, 4, blockedcrystal)
visualize(trsasnformedlocalregionnext, title="l", visualspace=euclidean(RealSpace, 2))

localeigspecnext = couriercorrelationRGthird[trsasnformedlocalfocknext, trsasnformedlocalfocknext] |> eigspech 
visualize(localeigspecnext)

#
physicalmodesRGthird = couriercorrelationRGthird |> getoutspace |> orderedmodes

refwannierfnthirdRG = locaclRGthird(center,couriercorrelationRGthird, physicalmodesRGthird, 0.006, 4)
localunitaryrefthirdRG = refwannierfnthirdRG[1]+refwannierfnthirdRG[2]

localunitarylistthirdRG = translatedlocalfockmaplist(center = center,localfockmap = localunitaryrefthirdRG,transceterlist=thirdlayerRGs)
wanniercourierlistthirdRG = translatedlocalfockmaplist(center = center,localfockmap = refwannierfnthirdRG[2],transceterlist=thirdlayerRGs)
localunitarythirdRG = sum(localunitarylistthirdRG)
wanniercouriersthirdRG = sum(wanniercourierlistthirdRG)

couriercorrelationRGRG = wanniercouriersthirdRG'*couriercorrelationRGthird[localunitarythirdRG |> getoutspace, localunitarythirdRG |> getoutspace]*wanniercouriersthirdRG

# inspecting local region spectrum
center = Point([0,0], scaledtriangular)
transformedphysicalmodesnext = couriercorrelationRGRG |> getoutspace |> orderedmodes

trsasnformedlocalregionnext, trsasnformedlocalfocknext = localregioninspection(center,transformedphysicalmodesnext, 8, blockedcrystal)
visualize(trsasnformedlocalregionnext, title="l", visualspace=euclidean(RealSpace, 2))

localeigspecnext = couriercorrelationRGRG[trsasnformedlocalfocknext, trsasnformedlocalfocknext] |> eigspech 
visualize(localeigspecnext, title="Spectrum after 2 RG steps for Dirac semi-metal")
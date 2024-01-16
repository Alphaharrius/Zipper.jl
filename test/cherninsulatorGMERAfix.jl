using Plotly, SmithNormalForm, LinearAlgebra, OrderedCollections, SparseArrays, Combinatorics
using Revise
using Zipper
using DataFrames
using SparseArrays

function sitefock(site::Offset; flavorcount::Integer = 1)::FockSpace{Offset}
    basis::Offset = site|>basispoint
    offset::Offset = site - basis
    return FockSpace((Mode(:offset => offset, :pos => basis, :flavor => f) for f in 1:flavorcount), reflected=site)
end

regionfock(region::Region; flavorcount::Integer = 1)::FockSpace{Region} = (
    fockspaceunion(sitefock(r, flavorcount=flavorcount) for r in region)|>FockSpace{Region})

function circularregionmodes(correlations::FockMap, origin::Offset, physicalmodes::Subset{Mode}, radius::Number)::Subset{Mode}
    currentspace::RealSpace = correlations.outspace |> getcrystal |> getspace |> orthospace
    physicalnorm = m -> lineartransform(currentspace, m |> getpos) |> norm
    return filter(p -> physicalnorm(p - origin) < radius, physicalmodes)
end

function circularregionmodesphysical(correlations::FockMap, origin::Offset,physicalmodes::Subset{Mode}, radius::Number)::Subset{Mode}
    currentspace::RealSpace = correlations.outspace |> getcrystal |> getspace |> euclidean
    physicalpos = m -> lineartransform(currentspace, m |> getpos) 
    return filter(p -> norm(physicalpos(p)-origin) < radius, physicalmodes)
end

function realphyspos(mode::Mode)
    currentspace::RealSpace = mode|> getspace |> euclidean
    return(lineartransform(currentspace, mode |> getpos))
end

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

function groupmodesbydistphysical(;
    region::Subset{Point{RealSpace}},
    regionfock::FockSpace,
    center::Point,
    samedistancethreshold::Int = 4)

    visualspace = region |> getspace |> euclidean
    physicalpos = m -> lineartransform(visualspace, m |> getpos) 
    distancewithmode = sort([(norm(physicalpos(mode)-center),mode) for mode in regionfock], by = first, rev = true)
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

function regioncorrelationsmap(correlations::FockMap, regionfock::FockSpace)::FockMap
    fouriermap::FockMap = fourier(correlations.inspace, regionfock) / (correlations.inspace |> subspacecount |> sqrt)
    return fouriermap
end


function givenewscale(scale::Scale, crystal::Crystal)
    realspace::RealSpace = crystal |> getspace
    snfinput::Matrix{Integer} = map(Integer, vcat(scale |> rep, crystal |> size |> diagm))
    U, S, Vd = snfinput |> dosnf
    snfinput = U[end - dimension(realspace) + 1:end, 1:dimension(realspace)]
    _, bS, bVd = snfinput |> dosnf
    newsizes = bS |> diag
    newrelativebasis = bVd * (S |> diag |> diagm) * Vd |> transpose

    subunitcell = Subset((transpose(Vd) * (r |> collect)) ∈ realspace for r in Iterators.product((0:(s-1) for s in S |> diag)...))
    newscale = Scale(newrelativebasis, realspace)
    newunitcell = Subset(newscale * (a + b) |> basispoint for (a, b) in Iterators.product(subunitcell, crystal |> getunitcell))
    return newscale
end

latticevec = [sqrt(3)/2 -1/2; 0. 1.]'
triangular = RealSpace([sqrt(3)/2 -1/2; 0. 1.]')
kspace = convert(MomentumSpace, triangular)

pa = [1/3, 2/3] ∈ triangular
pb = [2/3, 1/3] ∈ triangular
pc = (pa + pb) / 2
spatialsnappingcalibration((pa, pb, pc))

c6 = pointgrouptransform([cos(π/3) -sin(π/3); sin(π/3) cos(π/3)])

unitcell = Subset(pa, pb)
crystal = Crystal(unitcell, [24, 24])

reciprocalhashcalibration(crystal.sizes)

modes::Subset{Mode} = quantize(:pos, unitcell, 1)
m0, m1 = members(modes)

tₙ = -1 + 0im
tₕ = 0.1im

nearestneighbor = [
    (m0, m1) => tₙ,
    (m0, m1 |> setattr(:offset => [-1, 0] ∈ triangular)) => tₙ,
    (m0, m1 |> setattr(:offset => [0, 1] ∈ triangular)) => tₙ]

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
correlations = C
CC = CrystalFockMap(C)

crystalfock = correlations.outspace

scale = Scale([2 0; 0 2], crystalfock |> getcrystal |> getspace)
@time blockresult = blocking(:action => scale, :correlations => correlations, :crystal => (crystalfock |> getcrystal))

blockedcrystal::Crystal = blockresult[:crystal]
blockedcorrelations::FockMap = blockresult[:correlations]
blocker = blockresult[:transformer]

crystalpoints::Subset{Offset} = latticepoints(blockedcrystal)
samplepoints::Subset{Offset} = crystalpoints + c6^2 * crystalpoints + c6^4 * crystalpoints
blockedmodes::Subset{Mode} = quantize(:pos, blockedcrystal.unitcell, 1)
physicalmodes::Subset{Mode} = spanoffset(blockedmodes, samplepoints)

localmodes = circularregionmodes(blockedcorrelations, [1, 0] ∈ (blockedmodes |> getspace) , physicalmodes, 1)
localfock = FockSpace{Region}(localmodes)
visualize(getregion(localfock),visualspace=euclidean(RealSpace, 2))
visualize(regioncorrelations(blockedcorrelations, localfock) |> eigspec)

localcenter1::Offset = [1, 0] ∈ (blockedmodes |> getspace)
localmodes1 = circularregionmodes(blockedcorrelations, localcenter1 , physicalmodes, 1)
localfock1::FockSpace = FockSpace{Region}(localmodes1)
localregion1::Subset{Offset} = Subset(m |> getpos for m in localmodes1)
localregion1map = regioncorrelationsmap(blockedcorrelations, localfock1)
regioncorrelation1::FockMap = regioncorrelations(blockedcorrelations, localfock1)
modebydist1 = groupmodesbydist(region = localregion1,regionfock = localfock1,center = localcenter1)  
# group the eigenvectors to isometries with filled, empty ,courier
localiso1 = localisometries(blockedcorrelations, localfock1, selectionstrategy=modeselectionbycount(3))
# finding seeds for local frozen (6 modes at the center)
frozenseedingfocks1 = (FockSpace{Region}(circularregionmodes(blockedcorrelations, modewithdist[2]|> getpos, physicalmodes, 0.1)) for modewithdist in modebydist1[3])
frozenseeds1 = sum([reduce(+ ,FockMap(v,inspace = v |> getinspace |> orderedmodes |> setattr(:dumind => n) |> FockSpace,performpermute = false) for (n, (_ , v)) in enumerate(_findlocalspstates(statecorrelations = blockedcorrelations, regionfock = frozenseedingfock, 
spectrumextractpredicate = v -> true, symmetry = identitytransform(2)))) for frozenseedingfock in frozenseedingfocks1])
extendedfrozenseeds1 = idmap(localfock1, localfock1)[:,frozenseeds1 |> getoutspace] * frozenseeds1
wannierizedfrozens1 = _localwannierization(localiso1[:filled]+localiso1[:empty], extendedfrozenseeds1)

visualize(wannierizedfrozens1'*regioncorrelation1*wannierizedfrozens1 |> eigspec)

scaleforgmera = Scale([1 2; 2 1], blockedcrystal |> getspace)
visualize((scaleforgmera*blockedcrystal)|>getunitcell)

blockedbz = blockedcrystal|>brillouinzone
scaledbrillouinzone = scaleforgmera*blockedbz
scaledk = Subset(k|>basispoint for k in scaledbrillouinzone)
scaledk_dict = Dict()
for kpt in blockedbz
    key = basispoint(scaleforgmera*kpt)
    if haskey(scaledk_dict,key)
        push!(scaledk_dict[key],kpt)
    else
        scaledk_dict[key] = [kpt]
    end
end
scaledcrystal = scaleforgmera*blockedcrystal
scaledunitcell = scaledcrystal|> getunitcell
newscale = givenewscale(scaleforgmera,blockedcrystal)
unscaledunicell = inv(newscale)*scaledunitcell
newunitcellft = fourier(blockedcorrelations.inspace,regionfock(unscaledunicell))
fourier_subspaces = Dict(newunitcellft|>getoutspace|>crystalsubspaces)
scaled_fourier_list = []
for (scaled_k,old_ks) in scaledk_dict
    combined_space = fockspaceunion(fourier_subspaces[k] for k in old_ks)
    scaled_fourier_block = newunitcellft[combined_space,:]
    new_inspace = FockSpace(mode|>setattr(:pos=>scaleforgmera*getpos(mode),:offset=>scaled_k) for mode in scaled_fourier_block.inspace)
    scaled_fourier_block = FockMap(scaled_fourier_block,inspace = new_inspace,performpermute=false)
    push!(scaled_fourier_list,scaled_fourier_block)
end
scaled_fourier = directsum(scaled_fourier_list)
scaled_crystal_fock = FockSpace(scaled_fourier.inspace,reflected=scaledcrystal)
scaled_fourier = FockMap(scaled_fourier,inspace = scaled_crystal_fock,outspace = blockedcorrelations.inspace)
scaledcorrelations = scaled_fourier'*blockedcorrelations*scaled_fourier

scaledcrystalpoints::Subset{Offset} = latticepoints(scaledcrystal)
scaledsamplepoints::Subset{Offset} = scaledcrystalpoints + c6^2 * scaledcrystalpoints + c6^4 * scaledcrystalpoints
scaledmodes::Subset{Mode} = quantize(:pos, scaledcrystal.unitcell, 1)
scaledphysicalmodes::Subset{Mode} = spanoffset(scaledmodes, scaledsamplepoints)

localcenter = latticevec*[2,0]∈ (scaledphysicalmodes |> getspace) 
scaledlocalmodes = circularregionmodesphysical(scaledcorrelations, localcenter, scaledphysicalmodes, 2)
scaledlocalfock = FockSpace{Region}(scaledlocalmodes)
visualize(getregion(scaledlocalfock),visualspace=euclidean(RealSpace, 2))
visualize(regioncorrelations(scaledcorrelations, scaledlocalfock) |> eigspec)
localregion::Subset{Offset} = Subset(m |> getpos for m in scaledlocalmodes)
modebydist = groupmodesbydistphysical(region = localregion,regionfock = scaledlocalfock,center = localcenter)  

# group the eigenvectors to isometries with filled, empty ,courier
localiso = localisometries(scaledcorrelations, scaledlocalfock, selectionstrategy=modeselectionbycount(3))
# finding seeds for local frozen (6 modes at the center)
frozenseedingfocks = (FockSpace{Region}(circularregionmodesphysical(scaledcorrelations, realphyspos(modewithdist[2]), scaledphysicalmodes, 0.1)) for modewithdist in modebydist[3])
frozenseeds = sum([reduce(+ ,FockMap(v,inspace = v |> getinspace |> orderedmodes |> setattr(:dumind => n) |> FockSpace,performpermute = false) for (n, (_ , v)) in enumerate(_findlocalspstates(statecorrelations = scaledcorrelations, regionfock = frozenseedingfock, 
spectrumextractpredicate = v -> true, symmetry = identitytransform(2)))) for frozenseedingfock in frozenseedingfocks])
extendedfrozenseeds = idmap(scaledlocalfock, scaledlocalfock)[:,frozenseeds |> getoutspace] * frozenseeds
wannierizedfrozens = _localwannierization(localiso[:filled]+localiso[:empty], extendedfrozenseeds)

wanniercrystalisos = _crystalisometries(localisometry=wannierizedfrozens, crystalfock=scaledcorrelations.outspace, addinspacemomentuminfo = true)
wannierunitcell::Subset{Offset} = Subset((mode |> getattr(:pos))  for mode in frozenseeds.inspace |> orderedmodes)
visualize(wannierunitcell)

wanniercrystall::Crystal = Crystal(wannierunitcell, (scaledcorrelations.outspace |> getcrystal).sizes)

wannierisometry::FockMap = directsum(wanniercrystaliso for (k,wanniercrystaliso) in wanniercrystalisos)
inspace::FockSpace = FockSpace(wannierisometry.inspace, reflected=wanniercrystall)
outspace::FockSpace = FockSpace(wannierisometry.outspace, reflected=scaledcorrelations.outspace)
wannierizedfrozen = FockMap(wannierisometry, inspace=inspace, outspace=outspace, performpermute=false)

visualize((wannierizedfrozen'*scaledcorrelations*wannierizedfrozen)|>eigspec)

frozencorrelations = wannierizedfrozen'*scaledcorrelations*wannierizedfrozen
frozencrystal = getcrystal(frozencorrelation.inspace)
visualize(frozencrystal.unitcell)
frozencrystalpoints::Subset{Offset} = latticepoints(getcrystal(frozencorrelation.inspace))
frozensamplepoints::Subset{Offset} = frozencrystalpoints + c6^2 * frozencrystalpoints + c6^4 * frozencrystalpoints
frozenmodes::Subset{Mode} = quantize(:pos, frozencrystal.unitcell, 1)
frozenmodes|>modeattrs
frozenphysicalmodes::Subset{Mode} = spanoffset(frozenmodes, frozensamplepoints)

modeattrs(frozencorrelations.inspace)

localcenter = latticevec*[2,0]∈ (frozenphysicalmodes |> getspace) 
frozenlocalmodes = circularregionmodesphysical(frozencorrelations, localcenter, frozenphysicalmodes, 2)
frozenlocalfock = FockSpace{Region}(frozenlocalmodes)
visualize(getregion(frozenlocalfock),visualspace=euclidean(RealSpace, 2))
visualize(regioncorrelations(frozencorrelations, frozenlocalfock) |> eigspec)

(frozencorrelations.inspace) |> modeattrs

frozenlocalfock|> modeattrs
fourier(frozencorrelations.inspace, frozenlocalfock)

visualize(_regioncorrelations(frozencorrelations, frozenlocalfock) )

function _regioncorrelations(correlations::FockMap, regionfock::FockSpace)::FockMap
    fouriermap::FockMap = fourier(correlations.inspace, regionfock) / (correlations.inspace |> subspacecount |> sqrt)
    return fouriermap
end
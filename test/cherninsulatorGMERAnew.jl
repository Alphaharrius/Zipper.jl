using Plotly, SmithNormalForm, LinearAlgebra, OrderedCollections, SparseArrays, Combinatorics
using Revise
using Zipper
using DataFrames
using SparseArrays
# using AbstractAlgebra

# function Base.:*(scale::Scale, crystal::Crystal)::Crystal
#     realspace::RealSpace = crystal |> getspace
#     snfinput::Matrix{Integer} = map(Integer, vcat(scale |> rep, crystal |> size |> diagm))
#     U, S, Vd = snfinput |> dosnf
#     snfinput = U[end - dimension(realspace) + 1:end, 1:dimension(realspace)]
#     _, bS, bVd = snfinput |> dosnf
#     newsizes = bS |> diag
#     newrelativebasis = bVd * (S |> diag |> diagm) * Vd |> transpose

#     subunitcell = Subset((transpose(Vd) * (r |> collect)) ∈ realspace for r in Iterators.product((0:(s-1) for s in S |> diag)...))
#     newscale = Scale(newrelativebasis, realspace)
#     newunitcell = Subset(newscale * (a + b) |> basispoint for (a, b) in Iterators.product(subunitcell, crystal |> getunitcell))
#     return Crystal(newunitcell, newsizes)
# end

# function Zipper.fourier(crystal::Crystal, region::Region)
#     bz::Subset{Momentum} = crystal|>brillouinzone
#     barepoint::Offset = crystal|>getspace|>getorigin
#     barecrystal::Crystal = Crystal(barepoint|>Subset, crystal|>size)
#     outspace::CrystalFock = FockSpace((Mode(:offset => k, :pos => barepoint, :flavor => 1) for k in bz), reflected=barecrystal)
    
#     latticesites::Region = Subset(p - (p|>basispoint) for p in region) # Removing all sub-lattice degrees of freedom.
#     inspace::FockSpace{Region} = latticesites|>regionfock

#     momentummatrix::Matrix = hcat((k|>euclidean|>vec for k in bz)...)
#     offsetmatrix::Matrix = hcat((p|>euclidean|>vec for p in latticesites)...)

#     fouriermatrix::SparseMatrixCSC = exp.(-1im * momentummatrix' * offsetmatrix)
#     return FockMap(outspace, inspace, fouriermatrix)
# end

function sitefock(site::Offset; flavorcount::Integer = 1)::FockSpace{Offset}
    basis::Offset = site|>basispoint
    offset::Offset = site - basis
    return FockSpace((Mode(:offset => offset, :pos => basis, :flavor => f) for f in 1:flavorcount), reflected=site)
end

function tensorproduct(primary::FockSpace, secondary::FockSpace)
    return fockspaceunion(FockSpace(merge(pmode, smode) for smode in secondary) for pmode in primary)
end

Base.:merge(a::Mode, b::Mode)::Mode = Mode(Dict(getattrs(a)..., getattrs(b)...))

regionfock(region::Region; flavorcount::Integer = 1)::FockSpace{Region} = (
    fockspaceunion(sitefock(r, flavorcount=flavorcount) for r in region)|>FockSpace{Region})


function circularregionmodes(correlations::FockMap, origin::Offset, physicalmodes::Subset{Mode}, radius::Number)::Subset{Mode}
    currentspace::RealSpace = correlations.outspace |> getcrystal |> getspace |> orthospace
    physicalnorm = m -> lineartransform(currentspace, m |> getpos) |> norm
    return filter(p -> physicalnorm(p - origin) < radius, physicalmodes)
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

function _crystalisometries(; localisometry::FockMap, crystalfock::FockSpace{Crystal},
    addinspacemomentuminfo::Bool = false)::Dict{Momentum, FockMap}

    crystal::Crystal = getcrystal(crystalfock)
    # fouriermap::FockMap = fourier(crystalfock, localisometry.outspace) / (crystal |> vol |> sqrt) seems do not need the denominator
    fouriermap::FockMap = fourier(crystalfock, localisometry.outspace) 
    momentumfouriers::Base.Generator = rowsubmaps(fouriermap)
    bz::Subset{Momentum} = brillouinzone(crystal)


    function preprocesslocalisometry(k::Momentum)::FockMap
        if !addinspacemomentuminfo
            return localisometry
        end
        inspace::FockSpace = localisometry.inspace |> orderedmodes |> setattr(:koffset => k) |> FockSpace
        return FockMap(localisometry, inspace=inspace, performpermute=false)
    end

    return Dict(k => kfourier * preprocesslocalisometry(k) for (k, kfourier) in zip(bz, momentumfouriers))
end

function _crystalisometries(; localisometry::FockMap, crystalfock::FockSpace{Crystal},
    addinspacemomentuminfo::Bool = false)::Dict{Momentum, FockMap}

    crystal::Crystal = getcrystal(crystalfock)
    # fouriermap::FockMap = fourier(crystalfock, localisometry.outspace) / (crystal |> vol |> sqrt) seems do not need the denominator
    fouriermap::FockMap = fourier(crystalfock, localisometry.outspace) 
    momentumfouriers::Base.Generator = rowsubmaps(fouriermap)
    bz::Subset{Momentum} = brillouinzone(crystal)


    function preprocesslocalisometry(k::Momentum)::FockMap
        if !addinspacemomentuminfo
            return localisometry
        end
        inspace::FockSpace = localisometry.inspace |> orderedmodes |> setattr(:offset => k) |> FockSpace
        return FockMap(localisometry, inspace=inspace, performpermute=false)
    end

    return Dict(k => kfourier * preprocesslocalisometry(k) for (k, kfourier) in zip(bz, momentumfouriers))
end


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
@time blockresult = blocking(:action => scale, :correlations => correlations, :crystal => crystalfock |> getcrystal)

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
modebydist1 = groupmodesbydist(region = localregion1,regionfock = localfock1,center = localcenter1)  
# group the eigenvectors to isometries with filled, empty ,courier
localiso1 = localisometries(blockedcorrelations, localfock1, selectionstrategy=modeselectionbycount(3))
# finding seeds for local frozen (6 modes at the center)
frozenseedingfocks1 = (FockSpace{Region}(circularregionmodes(blockedcorrelations, modewithdist[2]|> getpos, physicalmodes, 0.1)) for modewithdist in modebydist1[3])
frozenseeds1 = sum([reduce(+ ,FockMap(v,inspace = v |> getinspace |> orderedmodes |> setattr(:dumind => n) |> FockSpace,performpermute = false) for (n, (_ , v)) in enumerate(_findlocalspstates(statecorrelations = blockedcorrelations, regionfock = frozenseedingfock, 
spectrumextractpredicate = v -> true, symmetry = identitytransform(2)))) for frozenseedingfock in frozenseedingfocks1])
extendedfrozenseeds1 = idmap(localfock1, localfock1)[:,frozenseeds1 |> getoutspace] * frozenseeds1
wannierizedfrozens1 = _localwannierization(localiso1[:filled]+localiso1[:empty], extendedfrozenseeds1)

scaleforgmera = Scale([1 2; 2 1], blockedcrystal |> getspace)
visualize((scaleforgmera*blockedcrystal)|>getunitcell)
scaleforgmera*blockedcorrelations.outspace

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
scaledcorrelation = scaled_fourier'*blockedcorrelations*scaled_fourier

# function getsphericalregion(; crystal::Crystal, radius::Real, metricspace::RealSpace)
#     generatingradius::Integer = ceil(Int, radius * 1.5) # Multiply by 1.5 to ensure all unitcell points fits.
#     generatinglength::Integer = generatingradius * 2
#     generatingcrystal::Crystal = Crystal(crystal|>getunitcell, [generatinglength, generatinglength])
#     crystalregion::Region = generatingcrystal|>sitepoints
#     centeredregion::Region = crystalregion - (crystalregion|>getcenter)
#     println(crystalregion)
#     println(centeredregion)
#     return Subset(point for point in centeredregion if norm(metricspace * point) <= radius)
# end

correlations

visualize(getregion(scaledlocalfock),visualspace=euclidean(RealSpace, 2))
visualize(regioncorrelations(scaledcorrelation, scaledlocalfock) |> eigspec)


blockedcorrelations.inspace |> modeattrs
regionfock(unscaledunicell) |> modeattrs

visualize(blockedbz, Subset(k|>basispoint for k in scaledbrillouinzone))

wannierizedfrozens1.inspace
#Construct crystal with this unitcell
localregion = localfock1|>getregion
unitcell = Subset(pt|>basispoint for pt in scaleforgmera*localregion)
gmeracrystal = Crystal(unitcell,[24,8])

crystalmap(blockedcorrelations.inspace,gmeracrystal)

crystalfock = blockedcorrelations.inspace
crystal = crystalfock |> getcrystal
scaledcrystal= gmeracrystal
kspace = scaledcrystal|>getspace
scale = (scaledcrystal|>getspace)*Scale(scaledcrystal|>getspace|>rep, euclidean(RealSpace,2))
unscaledblockedregion::Subset{Offset} = (scale |> inv) * scaledcrystal.unitcell
bz::Subset{Momentum} = crystal |> brillouinzone
basismodes::Subset{Mode} = crystalfock |> rep |> first
scaledbz::Subset{Momentum} = scaledcrystal |> brillouinzone

momentummappings::Base.Generator = (basispoint(scale * p) => p for p in bz)
mappingpartitions::Dict{Point, Vector{Point}} = foldl(momentummappings; init=Dict{Point, Vector{Point}}()) do d, (k, v)
        mergewith!(append!, d, LittleDict(k => [v]))
    end
ksubsets::Dict{Momentum, Subset{Mode}} = crystalfock |> crystalsubsets
scaledfock::FockSpace = ((ksubsets[k] for k in mappingpartitions[kscaled]) |> subsetunion |> FockSpace
    for kscaled in scaledbz) |> fockspaceunion
ksubsets

for kscaled in scaledbz
    print(kscaled)
end


(triangular |> rep) * (localregion|>getspace|>rep)
(triangular |> rep) *((scaleforgmera*localregion)|>getspace|>rep) 

# localcenter2::Offset = [0, 1] ∈ (blockedmodes |> getspace)
# localmodes2 = circularregionmodes(blockedcorrelations, localcenter2 , physicalmodes, 1)
# localfock2::FockSpace = FockSpace{Region}(localmodes2)
# localregion2::Subset{Offset} = Subset(m |> getpos for m in localmodes2)
# modebydist2 = groupmodesbydist(region = localregion2,regionfock = localfock2,center = localcenter2)  
# # group the eigenvectors to isometries with filled, empty ,courier
# localiso2 = localisometries(blockedcorrelations, localfock2, selectionstrategy=modeselectionbycount(3))
# # finding seeds for local frozen (6 modes at the center)
# frozenseedingfocks2 = (FockSpace{Region}(circularregionmodes(blockedcorrelations, modewithdist[2]|> getpos, physicalmodes, 0.1)) for modewithdist in modebydist2[3])
# frozenseeds2 = sum([reduce(+ ,FockMap(v,inspace = v |> getinspace |> orderedmodes |> setattr(:dumind => n) |> FockSpace,performpermute = false) for (n, (_ , v)) in enumerate(_findlocalspstates(statecorrelations = blockedcorrelations, regionfock = frozenseedingfock, 
# spectrumextractpredicate = v -> true, symmetry = identitytransform(2)))) for frozenseedingfock in frozenseedingfocks2])
# extendedfrozenseeds2 = idmap(localfock2, localfock2)[:,frozenseeds2 |> getoutspace] * frozenseeds2
# wannierizedfrozens2 = _localwannierization(localiso2[:filled]+localiso2[:empty], extendedfrozenseeds2)

# localcenter3::Offset = [-1, -1] ∈ (blockedmodes |> getspace)
# localmodes3 = circularregionmodes(blockedcorrelations, localcenter3 , physicalmodes, 1)
# localfock3::FockSpace = FockSpace{Region}(localmodes3)
# localregion3::Subset{Offset} = Subset(m |> getpos for m in localmodes3)
# modebydist3 = groupmodesbydist(region = localregion3,regionfock = localfock3,center = localcenter3)
# # group the eigenvectors to isometries with filled, empty ,courier
# localiso3 = localisometries(blockedcorrelations, localfock3, selectionstrategy=modeselectionbycount(3))
# # finding seeds for local frozen (6 modes at the center)
# frozenseedingfocks3 = (FockSpace{Region}(circularregionmodes(blockedcorrelations, modewithdist[2]|> getpos, physicalmodes, 0.1)) for modewithdist in modebydist3[3])
# frozenseeds3 = sum([reduce(+ ,FockMap(v,inspace = v |> getinspace |> orderedmodes |> setattr(:dumind => n) |> FockSpace,performpermute = false) for (n, (_ , v)) in enumerate(_findlocalspstates(statecorrelations = blockedcorrelations, regionfock = frozenseedingfock, 
# spectrumextractpredicate = v -> true, symmetry = identitytransform(2)))) for frozenseedingfock in frozenseedingfocks3])
# extendedfrozenseeds3 = idmap(localfock3, localfock3)[:,frozenseeds3 |> getoutspace] * frozenseeds3
# wannierizedfrozens3 = _localwannierization(localiso3[:filled]+localiso3[:empty], extendedfrozenseeds3)

# frozenseeds = frozenseeds1 + frozenseeds2 + frozenseeds3
# wannierizedfrozens = wannierizedfrozens1 + wannierizedfrozens2 + wannierizedfrozens3


# localfourier = Zipper.fourier(blockedcrystal, 
# FockSpace{Region}(Subset((mode)  for mode in frozenseeds.inspace |> orderedmodes))|> getregion)
# visualize(localfourier.inspace |> getregion)

# visualize(blockedcrystal |> getunitcell)
# fsr = FockSpace{Region}(frozenseeds1.inspace)  |> getregion

# visualize(Subset(pt + offset for offset in f.inspace |> getregion for pt in fsr), FockSpace{Region}(frozenseeds.inspace) |> getregion)

# localfourierFS = FockSpace(mode |> setattr(:r=> (mode |> getattr(:offset))) |> removeattr(:offset) for mode in (frozenseeds |> getinspace))

# fourierinspace = tensorproduct(localfourier|>getinspace,localfourierFS)
# fourieroutspace = tensorproduct(localfourier|>getoutspace,localfourierFS)

# 24*24*8


# 24*24
# (48*48)/576

# trialwannierizedfrozens = wannierizedfrozens1 + wannierizedfrozens2 + wannierizedfrozens3

# visualize(Subset(m |> getpos for m in (trialwannierizedfrozens.inspace |> getmodes)))


# localiso1::Dict{Symbol} = localisometries(correlations, localfock1, selectionstrategy=modeselectionbycount(3))
# localiso2::Dict{Symbol} = localisometries(correlations, localfock2, selectionstrategy=modeselectionbycount(3))
# localiso3::Dict{Symbol} = localisometries(correlations, localfock3, selectionstrategy=modeselectionbycount(3))

# blockedcrystal = blockedcorrelations.outspace |> getcrystal
# scaleforgmera = Scale([3 0; 0 3], blockedcrystal |> getspace)
# @time scaleresult = blocking(:action => scaleforgmera , :correlations => blockedcorrelations, :crystal => blockedcrystal)

# scaledcrystal = scaleresult[:crystal]
# scaledcorrelations = scaleresult[:correlations]
# scaler = scaleresult[:transformer]


A = map(Integer,vcat(scaleforgmera |> rep,blockedcrystal.sizes|>diagm))
U,S,Vd = dosnf(A)
U*S*Vd

(S*Vd)[1:2,1:2]
inv(U[3:end,1:2])*blockedcrystal.sizes|>diagm

u_tilde,s_tilde,v_tilde = dosnf(U[3:end,1:2])
v_tilde*(S*Vd)[1:2,1:2]
v_tilde*inv(U[3:end,1:2])*blockedcrystal.sizes|>diagm

# scaledcrystal = scaleforgmera*blockedcrystal
# visualize(scaledcrystal|>getunitcell)
realspace = blockedcrystal|>getspace
newrelativebasis = v_tilde * (S |> diag |> diagm) * Vd |> transpose
subunitcell = Subset((transpose(v_tilde) * (r |> collect)) ∈ realspace for r in Iterators.product((0:(s-1) for s in S |> diag)...))
newscale = Scale(newrelativebasis, realspace)
newunitcell = Subset(newscale * (a + b) |> basispoint for (a, b) in Iterators.product(subunitcell, blockedcrystal |> getunitcell))
visualize(newunitcell)

U, S, Vd = snfinput |> dosnf
    snfinput = U[end - dimension(realspace) + 1:end, 1:dimension(realspace)]
    _, bS, bVd = snfinput |> dosnf
    newsizes = bS |> diag
    newrelativebasis = bVd * (S |> diag |> diagm) * Vd |> transpose


scalingmap = scaleforgmera*blockedcorrelations.outspace
scaledcorrelations = scalingmap*blockedcorrelations*scalingmap'
scaledcrystal = scaledcorrelations.inspace |> getcrystal

refcrystalpoints::Subset{Offset} = latticepoints(scaledcrystal)
refsamplepoints::Subset{Offset} = refcrystalpoints + c6^2 * refcrystalpoints + c6^4 * refcrystalpoints
refblockedmodes::Subset{Mode} = quantize(:pos, scaledcrystal.unitcell, 1)
refphysicalmodes::Subset{Mode} = spanoffset(refblockedmodes, refsamplepoints)

# localmodes = circularregionmodes(scaledcorrelations, [0, 0] ∈ (refblockedmodes |> getspace) , refphysicalmodes, 0.3)
# localfock = FockSpace{Region}(localmodes)
# visualize(getregion(localfock),visualspace=euclidean(RealSpace, 2))
# visualize(regioncorrelations(scaledcorrelations, localfock) |> eigspec)

localcenter1::Offset = [1/3, 0] ∈ (refblockedmodes |> getspace)
localmodes1 = circularregionmodes(scaledcorrelations, localcenter1 , refphysicalmodes, 1/3)
localfock1::FockSpace = FockSpace{Region}(localmodes1)
localregion1::Subset{Offset} = Subset(m |> getpos for m in localmodes1)
modebydist1 = groupmodesbydist(region = localregion1,regionfock = localfock1,center = localcenter1) 
# group the eigenvectors to isometries with filled, empty ,courier
localiso1 = localisometries(scaledcorrelations, localfock1, selectionstrategy=modeselectionbycount(3))
# finding seeds for local frozen (6 modes at the center)
frozenseedingfocks1 = (FockSpace{Region}(circularregionmodes(scaledcorrelations, modewithdist[2]|> getpos, refphysicalmodes, 0.01)) for modewithdist in modebydist1[3])
frozenseeds1 = sum([reduce(+ ,FockMap(v,inspace = v |> getinspace |> orderedmodes |> setattr(:dumind => n) |> FockSpace,performpermute = false) for (n, (_ , v)) in enumerate(_findlocalspstates(statecorrelations = scaledcorrelations, regionfock = frozenseedingfock, 
spectrumextractpredicate = v -> true, symmetry = identitytransform(2)))) for frozenseedingfock in frozenseedingfocks1])
extendedfrozenseeds1 = idmap(localfock1, localfock1)[:,frozenseeds1 |> getoutspace] * frozenseeds1
wannierizedfrozens1 = _localwannierization(localiso1[:filled]+localiso1[:empty], extendedfrozenseeds1)

localcenter2::Offset = [0, 1/3] ∈ (refblockedmodes |> getspace)
localmodes2 = circularregionmodes(scaledcorrelations, localcenter2 , refphysicalmodes, 1/3)
localfock2::FockSpace = FockSpace{Region}(localmodes2)
localregion2::Subset{Offset} = Subset(m |> getpos for m in localmodes2)
modebydist2 = groupmodesbydist(region = localregion2,regionfock = localfock2,center = localcenter2)  
# group the eigenvectors to isometries with filled, empty ,courier
localiso2 = localisometries(scaledcorrelations, localfock2, selectionstrategy=modeselectionbycount(3))
# finding seeds for local frozen (6 modes at the center)
frozenseedingfocks2 = (FockSpace{Region}(circularregionmodes(scaledcorrelations, modewithdist[2]|> getpos, refphysicalmodes, 0.01)) for modewithdist in modebydist2[3])
frozenseeds2 = sum([reduce(+ ,FockMap(v,inspace = v |> getinspace |> orderedmodes |> setattr(:dumind => n) |> FockSpace,performpermute = false) for (n, (_ , v)) in enumerate(_findlocalspstates(statecorrelations = scaledcorrelations, regionfock = frozenseedingfock, 
spectrumextractpredicate = v -> true, symmetry = identitytransform(2)))) for frozenseedingfock in frozenseedingfocks2])
extendedfrozenseeds2 = idmap(localfock2, localfock2)[:,frozenseeds2 |> getoutspace] * frozenseeds2
wannierizedfrozens2 = _localwannierization(localiso2[:filled]+localiso2[:empty], extendedfrozenseeds2)

localcenter3::Offset = [-1/3, -1/3] ∈ (refblockedmodes |> getspace)
localmodes3 = circularregionmodes(scaledcorrelations, localcenter3 , refphysicalmodes, 1/3)
localfock3::FockSpace = FockSpace{Region}(localmodes3)
localregion3::Subset{Offset} = Subset(m |> getpos for m in localmodes3)
modebydist3 = groupmodesbydist(region = localregion3,regionfock = localfock3,center = localcenter3)
# group the eigenvectors to isometries with filled, empty ,courier
localiso3 = localisometries(scaledcorrelations, localfock3, selectionstrategy=modeselectionbycount(3))
# finding seeds for local frozen (6 modes at the center)
frozenseedingfocks3 = (FockSpace{Region}(circularregionmodes(scaledcorrelations, modewithdist[2]|> getpos, refphysicalmodes, 0.01)) for modewithdist in modebydist3[3])
frozenseeds3 = sum([reduce(+ ,FockMap(v,inspace = v |> getinspace |> orderedmodes |> setattr(:dumind => n) |> FockSpace,performpermute = false) for (n, (_ , v)) in enumerate(_findlocalspstates(statecorrelations = scaledcorrelations, regionfock = frozenseedingfock, 
spectrumextractpredicate = v -> true, symmetry = identitytransform(2)))) for frozenseedingfock in frozenseedingfocks3])
extendedfrozenseeds3 = idmap(localfock3, localfock3)[:,frozenseeds3 |> getoutspace] * frozenseeds3
wannierizedfrozens3 = _localwannierization(localiso3[:filled]+localiso3[:empty], extendedfrozenseeds3)

# visualize((wannierizedfrozens3'*regioncorrelations(scaledcorrelations, localfock3)*wannierizedfrozens3) |> eigspec)

frozenseeds = frozenseeds1 + frozenseeds2 + frozenseeds3
wannierizedfrozens = wannierizedfrozens1 + wannierizedfrozens2 + wannierizedfrozens3


localfourier = Zipper.fourier(scaledcrystal, 
FockSpace{Region}(Subset((mode)  for mode in frozenseeds.inspace |> orderedmodes))|> getregion)
visualize(localfourier.inspace |> getregion)

visualize(scaledcrystal |> getunitcell)
fsr = FockSpace{Region}(frozenseeds1.inspace)  |> getregion

visualize(Subset(pt + offset for offset in f.inspace |> getregion for pt in fsr), FockSpace{Region}(frozenseeds.inspace) |> getregion)

localfourierFS = FockSpace(mode |> setattr(:r=> (mode |> getattr(:offset))) |> removeattr(:offset) for mode in (frozenseeds |> getinspace))

fourierinspace = tensorproduct(localfourier|>getinspace,localfourierFS)
fourieroutspace = tensorproduct(localfourier|>getoutspace,localfourierFS)
fouriermtrix = kron(localfourier |> rep,idmap(frozenseeds.inspace)|> rep)
wannierizedfrozen = FockMap(fourieroutspace,fourierinspace, fouriermtrix)

fourierinspace|> modeattrs |> first
FockSpace(mode |> setattr(:offset=> (mode |> getattr(:r))) |> removeattr(:r) for mode in (fourierinspace))

72*8*8

visualize(latticepoints(scaledcorrelations.inspace|>getcrystal))

wanniercrystalisos = _crystalisometries(localisometry=wannierizedfrozens, crystalfock=scaledcorrelations.outspace, addinspacemomentuminfo = true)
wannierunitcell::Subset{Offset} = Subset((mode |> getattr(:pos))  for mode in frozenseeds.inspace |> orderedmodes)
visualize(wannierunitcell)
visualize(latticepoints(scaledcorrelations.outspace |> getcrystal))

wanniercrystall::Crystal = Crystal(wannierunitcell, (scaledcorrelations.outspace |> getcrystal).sizes)
# visualize(latticepoints(wanniercrystall))

wannierisometry::FockMap = directsum(wanniercrystaliso for (k,wanniercrystaliso) in wanniercrystalisos)
inspace::FockSpace = FockSpace(wannierisometry.inspace, reflected=wanniercrystall)
outspace::FockSpace = FockSpace(wannierisometry.outspace, reflected=scaledcorrelations.outspace)
wannierizedfrozen = FockMap(wannierisometry, inspace=inspace, outspace=outspace, performpermute=false)

frozencorrelation = (wannierizedfrozen'*scaledcorrelations*wannierizedfrozen)
visualize(frozencorrelation|> eigspec)


function _regioncorrelations(correlations::FockMap, regionfock::FockSpace)::FockMap
    fouriermap::FockMap = fourier(correlations.inspace, regionfock) / (correlations.inspace |> subspacecount |> sqrt)
    return fouriermap
end

frozencorrelation.inspace
(frozencorrelation.inspace |> getcrystal).unitcell
visualize((frozencorrelation.inspace |> getcrystal).unitcell)
trialcrystalpoints::Subset{Offset} = latticepoints(frozencorrelation.inspace |> getcrystal)
trialsamplepoints::Subset{Offset} = trialcrystalpoints + c6^2 * trialcrystalpoints + c6^4 * trialcrystalpoints
trialblockedmodes::Subset{Mode} = quantize(:pos, (frozencorrelation.inspace |> getcrystal).unitcell, 1)
trialphysicalmodes::Subset{Mode} = spanoffset(trialblockedmodes, trialsamplepoints)
trialmodes = circularregionmodes(frozencorrelation, [1/3+7,7] ∈ (trialblockedmodes |> getspace), trialphysicalmodes, 1/3)
trialfock = FockSpace{Region}(trialmodes)
modeattrs(trialfock)
modeattrs(frozencorrelation.outspace)
visualize(getregion(trialfock),visualspace=euclidean(RealSpace, 2))
Matrix(_regioncorrelations(frozencorrelation, trialfock).rep)
visualize(regioncorrelations(frozencorrelation, trialfock) |> eigspec)



refcrystalpoints::Subset{Offset} = latticepoints(scaledcrystal)
refsamplepoints::Subset{Offset} = refcrystalpoints + c6^2 * refcrystalpoints + c6^4 * refcrystalpoints
refblockedmodes::Subset{Mode} = quantize(:pos, scaledcrystal.unitcell, 1)
refphysicalmodes::Subset{Mode} = spanoffset(refblockedmodes, refsamplepoints)

localmodes = circularregionmodes(scaledcorrelations, [0, 0] ∈ (refblockedmodes |> getspace) , refphysicalmodes, 0.3)
localfock = FockSpace{Region}(localmodes)
visualize(getregion(localfock),visualspace=euclidean(RealSpace, 2))
visualize(regioncorrelations(scaledcorrelations, localfock) |> eigspec)

24*24*8
8*8*(8*3*3)
(8*3*3)
(24*3-6*3*3)*8*8
(24*3-6*3*3)
12*12*8

24*24*8
4*4*(8*6*6)
8*6*6
(96*3-30*3*3)*4*4
(96*3-30*3*3)
6*6*8

24*24*8
4*4*(8*6*6)
8*6*6
(96*3-24*3*3)*4*4
(96*3-24*3*3)
12*12*8


circularregionmodes(frozenseedingcenter2 , physicalmodes, 4)
getregion(frozenseedingfock3)

x = 1000000
a = 0
for m = 1:100
    # println((1-cos(2*pi*m*x))/m)
    a+=(1-cos((2*pi*m*x)/100))/m
end
a
log(x)

18*8*8
72*8*8
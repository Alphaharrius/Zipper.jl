include("../src/spaces.jl")
include("../src/geometries.jl")
include("../src/quantum.jl")
include("../src/physical.jl")
include("../src/plotting.jl")
include("../src/zer.jl")
include("../src/transformations.jl")

using PlotlyJS, LinearAlgebra, OrderedCollections, SparseArrays, ColorTypes, SmithNormalForm
using ..Spaces, ..Geometries, ..Quantum, ..Physical, ..Zer, ..Transformations
using ..Plotting

triangular = RealSpace([sqrt(3)/2 -1/2; 0. 1.]')
kspace = convert(MomentumSpace, triangular)

unitcell = union(Point([1/3, 2/3], triangular), Point([2/3, 1/3], triangular))
crystal = Crystal(unitcell, [32, 32])
modes::Subset{Mode} = quantize(:pos, unitcell, 1)
fullfock::FockSpace = crystalfock(modes, crystal)

scale = Scale([2 0; 0 2])
recipient = Recipient(fullfock, :crystal => crystal)
blockmap = scale * recipient
newcrystal = scale * crystal
newrecipient = Recipient(blockmap.outspace, :crystal => newcrystal)
newblockmap = scale * newrecipient

tₙ = ComplexF64(-1.)
m0, m1 = members(modes)
bonds::FockMap = bondmap([
    (m0, m1) => tₙ,
    (m0, setattr(m1, :offset => Point([-1, 0], triangular))) => tₙ,
    (m0, setattr(m1, :offset => Point([0, 1], triangular))) => tₙ])

𝐻::FockMap = hamiltonian(crystal, bonds)
visualize(𝐻, title="Hamiltonian", rowrange=1:64, colrange=1:64)
filled::FockMap = groundstates(𝐻)

𝐶::FockMap = idmap(filled.outspace, filled.outspace) - filled * filled'
visualize(𝐶, title="Correlation", rowrange=1:64, colrange=1:64)
blockedcorrelation = blockmap * 𝐶 * blockmap'
visualize(blockedcorrelation, title="Correlation", rowrange=1:64, colrange=1:64)
firstfock = FockSpace(first(rep(blockedcorrelation.outspace)))
onepart = restrict(blockedcorrelation, firstfock, firstfock)
visualize(onepart)

newblockedcorrelation = newblockmap * blockedcorrelation * newblockmap'
visualize(newblockedcorrelation, rowrange=1:128, colrange=1:128)

newcrystal::Crystal = scale * crystal
blockedregion::Subset{Point} = inv(scale) * newcrystal.unitcell
[p for p in blockedregion]
BZ::Subset{Point} = brillouinzone(crystal)
fockspace::FockSpace = sparsefock(modes, BZ)
basismodes::Subset{Mode} = fockspace |> rep |> first
newBZ::Subset{Point} = brillouinzone(newcrystal)
# Generate the Dict which keys each fockspace partition by its momentum.
momentumtopartition::Dict{Point, Subset{Mode}} = Dict(getattr(first(part), :offset) => part for part in rep(fockspace))
momentummappings::Vector{Pair{Point, Point}} = [basispoint(scale * p) => p for p in BZ]
mappingpartitions::Dict{Point, Vector{Point}} = foldl(momentummappings; init=Dict{Point,Vector{Point}}()) do d,(k,v)
    mergewith!(append!, d, LittleDict(k=>[v]))
end
blockedcrystalpartitions::Subset{Subset{Mode}} = Subset([union([momentumtopartition[k] for k in mappingpartitions[scaled_k]]...) for scaled_k in newBZ])
blockedcrystalordering::Dict{Mode, Integer} = Dict(mode => index for (index, mode) in enumerate(flatten(blockedcrystalpartitions)))
blockedcrystalfock::FockSpace = FockSpace(blockedcrystalpartitions, blockedcrystalordering)
[m.attrs for m in [orderedmodes(blockedcrystalfock)...][1:33]]

blockedoffsets::Subset{Point} = Subset([pbc(crystal, p) |> latticeoff for p in blockedregion])
blockedfock::FockSpace = FockSpace(spanoffset(basismodes, blockedoffsets))
[m.attrs for m in orderedmodes(blockedfock)]

restrictedfourier::FockMap = fourier(BZ, blockedfock)
[m.attrs for m in [orderedmodes(restrictedfourier.outspace)...][1:33]]
visualize(restrictedfourier, rowrange=1:40)

permutedmap::FockMap = Quantum.permute(restrictedfourier, blockedcrystalfock, restrictedfourier.inspace) / sqrt(vol(crystal) / vol(newcrystal))
visualize(permutedmap, rowrange=1:40)
orderingrule(blockedcrystalfock, restrictedfourier.outspace)
function repack_fourierblocks(sourcemap::FockMap, scaled_k::Point, partition::Subset{Mode})::FockMap
    mappart::FockMap = rows(sourcemap, FockSpace(partition))
    inmodes::Subset{Mode} = Subset([
        setattr(m, :groups => ModeGroup(transformed, "scaled"), :index => i, :offset => scaled_k, :pos => convert(Point, m))
        for (i, m) in enumerate(orderedmodes(mappart.inspace))])
    return FockMap(mappart.outspace, FockSpace(inmodes), rep(mappart))
end
mapblocks::Vector{FockMap} = [repack_fourierblocks(permutedmap, scaled_k, partition) for (scaled_k, partition) in Iterators.zip(newBZ, rep(blockedcrystalfock))]
blockmap::FockMap = focksum(mapblocks)
visualize(blockmap)

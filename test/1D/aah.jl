include("../../src/spaces.jl")
include("../../src/quantum.jl")
include("../../src/physical.jl")
include("../../src/plotting.jl")
include("../../src/geometries.jl")
include("../../src/transformations.jl")
include("../../src/zer.jl")

using PlotlyJS
using ..Spaces, ..Quantum, ..Physical, ...Plotting, ...Geometries, ...Transformations, ...Zer

𝑅₁ = euclidean(RealSpace, 1)
unitcell::Subset{Point} = Subset(Point([1/2], 𝑅₁))
quantized::Subset{Mode} = quantize(:pos, unitcell, 1)
mode::Mode = first(quantized)

# Generate all modes spanning the space.
𝑁::Integer = 64 # Number of lattice sites.
crystal::Crystal = Crystal(unitcell, [𝑁])
modes::Subset{Mode} = spanoffset(Subset(mode), latticepoints(crystal))
[convert(Point, m) for m in modes]

# ======================================================================================================
# Generate the AAH hamiltonian.
𝑉::Float64 = 0
𝑡::Float64 = 1
α::Float64 = (√5 + 1) / 2
onsite_bonds::Dict{Tuple{Mode, Mode}, ComplexF64} = Dict()
foreach(tup -> onsite_bonds[(last(tup), last(tup))] = 𝑉 * cos(2 * π * α * first(tup)), enumerate(modes)) # Fill the diagonal terms.
nn_bonds::Dict{Tuple{Mode, Mode}, ComplexF64} = Dict()
foreach(tup -> nn_bonds[(last(tup), modes[first(tup) + 1])] = -𝑡, enumerate([modes...][1:𝑁-1])) # Fill the nearest neighbor terms.
𝐹ᵧ::FockSpace = FockSpace(modes)
𝐻ₙₙ::FockMap = FockMap(𝐹ᵧ, 𝐹ᵧ, nn_bonds)
𝐻::FockMap = FockMap(𝐹ᵧ, 𝐹ᵧ, onsite_bonds) + 𝐻ₙₙ + 𝐻ₙₙ'
# ======================================================================================================

# ===========================================================================================
# Compute the ground state correlations.
eigenvalues::Vector{Pair{Mode, Float64}} = eigvalsh(𝐻)
plot(scatter(y=map(p -> p.second, eigenvalues), mode="markers"))

eigenvectors::FockMap = eigvecsh(𝐻)
visualize(eigenvectors)

groundmodes::Subset{Mode} = Subset(map(p -> p.first, filter(p -> p.second < 0, eigenvalues)))
groundstates::FockMap = columns(eigenvectors, FockSpace(groundmodes))

𝐶::FockMap = groundstates * groundstates'
# ===========================================================================================

visualize(𝐻)
visualize(𝐶)

scale::Scale = Scale([2.0][:,:])
newcrystal::Crystal = scale * crystal
newbasismodes::Subset{Mode} = quantize(:pos, newcrystal.unitcell, 1)
newmodes::Subset{Mode} = spanoffset(newbasismodes, latticepoints(newcrystal))

modemapping::Dict{Point, Mode} = Dict(convert(Point, m) => m for m in newmodes)
allignment::Subset{Mode} = Subset(modemapping[scale * convert(Point, m)] for m in modes)
newlattfock::FockSpace = FockSpace(newmodes)
temp::FockMap = Quantum.permute(idmap(newlattfock, newlattfock), inspace=FockSpace(allignment))
permutation::FockMap = FockMap(temp.outspace, FockSpace(modes), rep(temp))

blockedcorrelations::FockMap = permutation * 𝐶 * permutation'

restrictedregion::Subset{Mode} = filter(circularfilter(origin(euclidean(RealSpace, 1)), 4.0), newmodes)
restrictedfock::FockSpace = FockSpace(restrictedregion)
restrictedcorrelations::FockMap = restrict(blockedcorrelations, restrictedfock, restrictedfock)
visualize(restrictedcorrelations)
rcvals, rcunitary = eigh(restrictedcorrelations)
plot(scatter(y=map(p -> p.second, rcvals), mode="markers"), Layout(title="Local Distiller 𝑉 = $(𝑉)"))

distillers::Vector{FockMap} = []
for offset in [latticepoints(newcrystal)...][3:end-1]
    distillregion::Subset{Mode} = restrictedregion + offset
    distillfock::FockSpace = FockSpace(distillregion)
    distillcorrelations::FockMap = restrict(blockedcorrelations, distillfock, distillfock)
    frozenfocks = frozenselectionbycount(1)(distillcorrelations)
    localdistiller = frozenfocks[:empty] * frozenfocks[:empty]' - frozenfocks[:filled] * frozenfocks[:filled]'
    push!(distillers, localdistiller)
end
globaldistiller::FockMap = focksum(distillers)
visualize(globaldistiller)
gdvals, globalunitary = eigh(globaldistiller)
plot(scatter(y=map(p -> p.second, gdvals), mode="markers"), Layout(title="Global Distiller 𝑉 = $(𝑉)"))
globalcourier::Subset{Mode} = Subset(map(p -> p.first, filter(p -> -0.12 < p.second < 0.12, gdvals)))
globalcourierisometry::FockMap = columns(globalunitary, FockSpace(globalcourier))
visualize(globalcourierisometry)
mode = orderedmodes(globalcourierisometry.inspace)[21]
modespec = columnspec(columns(globalunitary, FockSpace(Subset(mode))))
visualize(modespec)

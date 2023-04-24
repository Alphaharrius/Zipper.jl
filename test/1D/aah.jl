include("../../src/spaces.jl")
include("../../src/quantum.jl")
include("../../src/physical.jl")

using PlotlyJS
using ..Spaces, ..Quantum, ..Physical

𝑅₁ = euclidean(RealSpace, 1)
unitcell::Subset{Point} = Subset([Point([1/2], 𝑅₁)])
quantized::Subset{Mode} = quantize("physical", :pos, unitcell, 1)
mode::Mode = first(quantized)

# Generate all modes spanning the space.
𝑁::Integer = 64 # Number of lattice sites.
modes::Vector{Mode} = [setattr(mode, :offset => Point([n], 𝑅₁)) for n in 1:𝑁]

# Generate the AAH hamiltonian
𝑉::Float64 = 0
𝑡::Float64 = 1
α::Float64 = (√5 + 1) / 2
onsite_bonds::Dict{Tuple{Mode, Mode}, ComplexF64} = Dict()
foreach(tup -> onsite_bonds[(last(tup), last(tup))] = 𝑉 * cos(2 * π * α * first(tup)), enumerate(modes)) # Fill the diagonal terms.
nn_bonds::Dict{Tuple{Mode, Mode}, ComplexF64} = Dict()
foreach(tup -> nn_bonds[(last(tup), modes[first(tup) + 1])] = -𝑡, enumerate(modes[1:𝑁-1])) # Fill the nearest neighbor terms.
𝐹ᵧ::FockSpace = FockSpace(Subset(modes))
𝐻ₙₙ::FockMap = FockMap(𝐹ᵧ, 𝐹ᵧ, nn_bonds)
𝐻::FockMap = FockMap(𝐹ᵧ, 𝐹ᵧ, onsite_bonds) + 𝐻ₙₙ + 𝐻ₙₙ'

eigenvalues::Vector{Pair{Mode, Float64}} = eigvalsh(𝐻)
plot(scatter(y=map(p -> p.second, eigenvalues)))

eigenvectors::FockMap = eigvecsh(𝐻)
plot(heatmap(z=real(rep(eigenvectors))))

𝐶::FockMap = ground_state_correlation(𝐻)
plot(heatmap(z=real(rep(𝐶))))

block_fockspace::FockSpace = FockSpace(Subset(modes[1:32]))
𝐶ᵣ::FockMap = rows(columns(𝐶, block_fockspace), block_fockspace)
plot(heatmap(z=real(rep(𝐶ᵣ))))

restricted_corrvals::Vector{Pair{Mode, Float64}} = eigvalsh(𝐶ᵣ)
plot(scatter(y=map(p -> p.second, restricted_corrvals), mode="markers"), Layout(title="𝑉 = $(𝑉)"))

include("../src/spaces.jl")
include("../src/geometries.jl")
include("../src/quantum.jl")
include("../src/physical.jl")
include("../src/plotting.jl")

using LinearAlgebra, PlotlyJS, OrderedCollections, SparseArrays, ColorTypes 
using ..Spaces, ..Geometries, ..Quantum, ..Physical, ..Plotting

triangular = RealSpace([sqrt(3)/2 -1/2; 0. 1.]')

k_space = convert(MomentumSpace, triangular)

unitcell = union(Point([1/3, 2/3], triangular), Point([2/3, 1/3], triangular))
crystal = Crystal(unitcell, [32, 32])
modes::Subset{Mode} = quantize("physical", :pos, unitcell, 1)
m0, m1 = members(modes)
tₙ = ComplexF64(-1.)
bonds::FockMap = bondmap([
    (m0, m1) => tₙ,
    (m0, setattr(m1, :offset => Point([-1, 0], triangular))) => tₙ,
    (m0, setattr(m1, :offset => Point([0, 1], triangular))) => tₙ])

𝐻::FockMap = @time hamiltonian(crystal, bonds)
filled::FockMap = @time filledstates(𝐻)
𝐶::FockMap = filled * filled'
plot(heatmap(z=real(rep(𝐶))))
cvs = eigvalsh(𝐶)
plot(scatter(y=map(p -> p.second, cvs), mode="markers"))
𝔘::FockMap = eigvecsh(𝐶)
plot(heatmap(z=real(rep(𝔘))))

small_crystal = Crystal(unitcell, [2, 2])
restricted_fockspace = FockSpace(Subset([setattr(m, :offset => p) for p in latticepoints(small_crystal) for m in modes]))

𝐹::FockMap = @time fourier(brillouin_zone(crystal), restricted_fockspace) / sqrt(vol(crystal))
𝐶ᵣ::FockMap = 𝐹' * 𝐶 * 𝐹
plot(heatmap(z=real(rep(𝐶ᵣ))))
𝑈ᵣ::FockMap = eigvecsh(𝐶ᵣ)
plot(heatmap(z=imag(rep(𝑈ᵣ))))
emode::Mode = orderedmodes(𝑈ᵣ.inspace)[1]
moderep::FockMap = columns(𝑈ᵣ, FockSpace(Subset([emode])))
values = columnspec(moderep)

function vv(title::String, spectrum::Vector{Pair{Mode, T}}) where {T <: Number}
    function generatetrace(element::Pair{Mode, T}) where {T <: Number}
        offset::Point = getattr(element.first, :offset)
        pos::Point = linear_transform(spaceof(offset), getattr(element.first, :pos))
        point::Point = offset + pos
        position::Vector{Float64} = Spaces.pos(linear_transform(euclidean(RealSpace, dimension(point)), point))
        padded_position::Vector{Float64} = vcat(position, zeros(Float64, 3 - length(position)))
        value::ComplexF64 = ComplexF64(element.second)
        @show value
        # Generate the RGB color based on the complex value.
        hue = angle(value) / 2π
        rgb::RGB{Float64} = convert(RGB{Float64}, HSV(hue, 1, 1))
        colortext::String = "rgba($(rgb.r * 255), $(rgb.g * 255), $(rgb.b * 255), 0.3)"
        @show colortext
        return scatter3d(x=[padded_position[1]], y=[padded_position[2]], z=[padded_position[3]], mode="marker", marker=attr(
            size=abs(value) * 50,
            color=colortext))
    end

    layout::Layout = Layout(title=title, scene=attr(aspectmode="data"))
    plot([generatetrace(element) for element in spectrum], layout)
end

function vvv(title::String, spectrum::Vector{Pair{Mode, T}}) where {T <: Number}
    𝑁::Int64 = length(spectrum)
    ∑𝑝::Vector{Point} = [getattr(pair.first, :offset) + getattr(pair.first, :pos) for pair in spectrum]
    𝑀ₚ::Matrix{Float64} = hcat(map(𝑝 -> pos(linear_transform(euclidean(RealSpace, dimension(𝑝)), 𝑝)), ∑𝑝)...)
    markerpositions::Matrix{Float64} = zeros(3, 𝑁)
    copyto!(view(markerpositions, 1:size(𝑀ₚ, 1), :), 𝑀ₚ)
    sizes::Vector{Float64} = [abs(pair.second) for pair in spectrum]
    markersizes::Vector{Float64} = sizes / norm(sizes) * 80
    colors::Vector{RGB{Float64}} = [convert(RGB{Float64}, HSV(angle(pair.second) / 2π, 1, 1)) for pair in spectrum]
    markercolors::Vector{Tuple{Float32, Float32, Float32}} = map(c -> Tuple([c.r, c.g, c.b] * 255), colors)
    @show [pair.second for pair in spectrum]
    trace = scatter3d(x=markerpositions[1, :], y=markerpositions[2, :], z=markerpositions[3, :], mode="markers", marker=attr(
        size=markersizes,
        color=markercolors))
    layout::Layout = Layout(title=title, scene=attr(aspectmode="data"))
    plot([trace], layout)
end

vvv("Mode", values)

rng = -1:0.1:1
tspec::Matrix{Float64} = zeros(Float64, length(rng), length(rng))
bspec::Matrix{Float64} = zeros(Float64, length(rng), length(rng))
for (j, k) in enumerate(rng)
    for (i, h) in enumerate(rng)
        fkk::FockMap = fourier(Subset([Point([h, k], k_space)]), Subset(orderedmodes(bm.outspace)))
        bloch::FockMap = fkk * bm * fkk'
        spec = eigvalsh(bloch)
        tspec[j, i] = real(spec[1].second)
        bspec[j, i] = real(spec[2].second)
    end
end

# spec = [hcat([eigvalsh(bloch(bonds, Point([h, k], k_space), crystal), :offset => Point([h, k], k_space)) for h in -1:0.1:1]...) for k in -1:0.1:1]

# top = map(p -> linear_transform(euclidean(MomentumSpace, 2), getattr(p.first, :offset)) => p.second, hcat([v[1, :] for v in spec]...))
# bottom = map(p -> linear_transform(euclidean(MomentumSpace, 2), getattr(p.first, :offset)) => p.second, hcat([v[2, :] for v in spec]...))

# topvals = hcat([map(p -> [pos(p.first)..., p.second], top)...]...)
# botvals = hcat([map(p -> [pos(p.first)..., p.second], bottom)...]...)

plot([surface(z=tspec), surface(z=bspec)])



# visualize_region("Honeycomb lattice", real_zone, euclidean(RealSpace, 2))
# visualize_region("Honeycomb lattice", reciprocal_zone, euclidean(RealSpace, 2))

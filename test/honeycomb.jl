include("../src/spaces.jl")
include("../src/geometries.jl")
include("../src/quantum.jl")
include("../src/physical.jl")
include("../src/plotting.jl")

using LinearAlgebra, PlotlyJS
using ..Spaces, ..Geometries, ..Quantum, ..Physical, ..Plotting

triangular = RealSpace([sqrt(3)/2 -1/2; 0. 1.]')

k_space = convert(MomentumSpace, triangular)

unitcell = union(Point([1/3, 2/3], triangular), Point([2/3, 1/3], triangular))
crystal = Crystal(unitcell, [8, 8])
modes::Subset{Mode} = quantize("physical", :pos, unitcell, 1)
m0, m1 = members(modes)

tₙ = ComplexF64(-1.)
bm = bondmap([
    (m0, m1) => tₙ,
    (m0, setattr(m1, :offset => Point([-1, 0], triangular))) => tₙ,
    (m0, setattr(m1, :offset => Point([0, 1], triangular))) => tₙ])

bzone::Subset{Point} = brillouin_zone(crystal)

struct Obj
    x
    y
end

obj_arr = [Obj(1,2), Obj(3,4), Obj(5,6)]

get_xy(obj) = (obj.x, obj.y)

xs, ys = map(get_xy, obj_arr) |> Base.transpose |> Base.vec

function directsum(elements::Vector{FockMap})::FockMap
    outparts = Iterators.map(el -> orderedmodes(el.outspace), elements)
    inparts = Iterators.map(el -> orderedmodes(el.inspace), elements)
    lengths::Vector{Integer} = [sum(part -> length(part), outparts), sum(part -> length(part), inparts)]
    outmodes::Vector{Mode} = vcat([outparts...]...)
    inmodes::Vector{Mode} = vcat([inparts...]...)
    spmat::SparseMatrixCSC{ComplexF64, Int64} = spzeros(lengths...)
    outordering::Dict{Mode, Integer} = Dict(Iterators.map(tup -> first(tup) => last(tup), enumerate(outmodes)))
    inordering::Dict{Mode, Integer} = Dict(Iterators.map(tup -> first(tup) => last(tup), enumerate(inmodes)))
    for (spidx, fockpart) in enumerate(fockparts)
        rowslice = outordering[first(fockpart)]:outordering[last(fockpart)]
        colslice = inordering[first(fockpart)]:inordering[last(fockpart)]
        spmat[rowslice, colslice] = rep(elements[spidx])
    end
    outspace::FockSpace = FockSpace(Subset(Set(Iterators.map(part -> Subset(part), outparts))))
    inspace::FockSpace = FockSpace(Subset(Set(Iterators.map(part -> Subset(part), inparts))))
    return FockMap(outspace, inspace, spmat)
end

ks = [bzone...][1:3]

function make_bloch(k::Point)::FockMap
    fmap::FockMap = fourier(Subset([k]), Subset(orderedmodes(bm.outspace)))
    return fmap * bm * dagger(fmap)
end

fockmaps::Vector{FockMap} = [make_bloch(k) for k in ks]

outparts = Iterators.map(el -> orderedmodes(el.outspace), fockmaps)
vcat([outparts...]...)

directsum(fockmaps)

rng = -1:0.1:1
tspec::Matrix{Float64} = zeros(Float64, length(rng), length(rng))
bspec::Matrix{Float64} = zeros(Float64, length(rng), length(rng))
for (j, k) in enumerate(rng)
    for (i, h) in enumerate(rng)
        fkk::FockMap = fourier(Subset([Point([h, k], k_space)]), Subset(orderedmodes(bm.outspace)))
        bloch::FockMap = fkk * bm * dagger(fkk)
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

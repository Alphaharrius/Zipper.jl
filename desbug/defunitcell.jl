using Zipper
using LinearAlgebra,Plots

function basispoint(point::Point)::Point
    rationalized::Vector = [hashablereal(v) for v in point |> vec]
    return Point([mod(v |> numerator, v |> denominator) / denominator(v) for v in rationalized], point |> getspace)
end

triangular = RealSpace([1 0; 1/2 sqrt(3)/2]')
kspace = convert(MomentumSpace, triangular)

# pa = [1/3, 0] ∈ triangular
# pb = [0, 1/3] ∈ triangular
# pc = [0, 2/3] ∈ triangular
# pd = [1/3, 2/3] ∈ triangular
# pe = [2/3, 1/3] ∈ triangular
# pf = [2/3, 0] ∈ triangular

pa = [1/3, 0] ∈ triangular
pb = [0, 1/3] ∈ triangular
pc = [-1/3, 1/3] ∈ triangular
pd = [-1/3, 0] ∈ triangular
pe = [0, -1/3] ∈ triangular
pf = [1/3, -1/3] ∈ triangular

pg = (pa + pb + pc + pd + pe + pf) / 6
spatialsnappingcalibration((pa, pb, pc, pd, pe, pf, pg))

unitcell = Subset(pa, pb, pc, pd, pe, pf)
crystal = Crystal(unitcell, [systemsize, systemsize])

reciprocalhashcalibration(crystal.sizes)
    
modes::Subset{Mode} = quantize(unitcell, 1)|>orderedmodes
pc|>basispoint

[m|>getattr(:r) for m in modes]
m0, m1, m2, m3, m4, m5 = members(modes)

onsite = [
    (m0, m0) => t_a,
    (m1, m1) => t_b,
    (m2, m2) => t_a,
    (m3, m3) => t_b,
    (m4, m4) => t_a,
    (m5, m5) => t_b
]
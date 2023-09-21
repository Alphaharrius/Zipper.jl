include("../src/spaces.jl")
include("../src/geometries.jl")
include("../src/transformations.jl")
include("../src/quantum.jl")
include("../src/plotting.jl")

using LinearAlgebra, Combinatorics
using ..Spaces, ..Geometries, ..Transformations, ..Plotting

function computeeigenfunctions(pointgroup::AffineTransform; dimension::Integer = pointgroup |> dimension, functionorderrange::UnitRange = 1:3)::Base.Iterators.Flatten
    function computeeigenfunctionsatorder(functionorder::Integer)
        matrixatorder::Matrix = pointgrouprepresentation(pointgroup.transformmatrix; rank=functionorder)
        eigenvalues, eigenvectors = matrixatorder |> eigen
        basisfunctions = (eigenvalue => BasisFunction(eigenvectors[:, n], dimension, functionorder) |> normalize for (n, eigenvalue) in eigenvalues |> enumerate)
        return Iterators.filter(p -> !(p.second |> iszero), basisfunctions)
    end

    return (p for functionorder in functionorderrange for p in functionorder |> computeeigenfunctionsatorder)
end

function computecharactertable(pointgroup::AffineTransform, eigenfunctions::Base.Iterators.Flatten, eigenvaluehashdenominator::Integer)::Matrix
    classfunctions = Dict(hashablecomplex(eigenvalue, eigenvaluehashdenominator) => (eigenvalue, eigenfunction) for (eigenvalue, eigenfunction) in eigenfunctions) |> values
    sortedclassfunctions = sort(classfunctions |> collect, by=p -> p |> first |> angle |> abs)
    transformelements::Vector{AffineTransform} = pointgroup |> pointgroupelements
    computecharacterrow(classfunction::BasisFunction) = [(element * classfunction) |> relativephase(classfunction) for element in transformelements] |> transpose
    return vcat((classfunction |> computecharacterrow for (_, classfunction) in sortedclassfunctions)...)
end

struct QuantumAffine <: Transformation{Matrix}
    affine::AffineTransform
    eigenfunctions::Dict{Tuple, BasisFunction}
    charactertable::Matrix
    eigenvaluehashdenominator::Integer
end

function functionclass(affine::AffineTransform, charactertable::Matrix, basis::BasisFunction)::Integer
    phasevector::Vector = [normalize(element * basis) |> relativephase(basis |> normalize) for element in affine |> pointgroupelements]
    classvector::Vector = map(v -> isapprox(v, 0, atol=1e-6) ? 0 : v, charactertable * phasevector)
    return findfirst(abs.(classvector) .> 0)
end

functionsignature(basis::BasisFunction, class::Integer)::Tuple = (basis |> dimension, basis.rank, class)

function QuantumAffine(pointgroup::AffineTransform)
    eigenfunctions::Base.Iterators.Flatten = pointgroup |> computeeigenfunctions
    eigenvaluehashdenominator::Integer = ((v for (v, _) in eigenfunctions) |> findcomplexdenominator).denominator
    charactertable::Matrix = computecharactertable(pointgroup, eigenfunctions, eigenvaluehashdenominator)
    eigenfunctiontable::Dict{Tuple, BasisFunction} = Dict(functionsignature(f, functionclass(pointgroup, charactertable, f)) => f for (_, f) in eigenfunctions)
    return QuantumAffine(pointgroup, eigenfunctiontable, charactertable, eigenvaluehashdenominator)
end

Base.convert(::Type{AffineTransform}, source::QuantumAffine)::AffineTransform = source.affine

Transformations.recenter(quantumaffine::QuantumAffine, center::Position)::QuantumAffine = QuantumAffine(
    quantumaffine.affine |> recenter(center), quantumaffine.eigenfunctions, quantumaffine.charactertable, quantumaffine.eigenvaluehashdenominator)
Transformations.affinematrix(quantumaffine::QuantumAffine)::Matrix = quantumaffine.affine |> affinematrix
Spaces.getspace(quantumaffine::QuantumAffine)::Space = quantumaffine.affine |> getspace
Base.:convert(::Type{Matrix}, source::QuantumAffine)::Matrix = source |> affinematrix
Base.:(==)(a::QuantumAffine, b::QuantumAffine)::Bool = a.affine == b.affine
Base.:*(space::RealSpace, quantumaffine::QuantumAffine)::QuantumAffine = QuantumAffine(
    space * quantumaffine.affine, quantumaffine.eigenfunctions, quantumaffine.charactertable, quantumaffine.eigenvaluehashdenominator)
Base.:*(a::QuantumAffine, b::QuantumAffine)::QuantumAffine = QuantumAffine(a.affine * b.affine)
Base.:*(quantumaffine::QuantumAffine, space::RealSpace)::QuantumAffine = quantumaffine.affine * space
Base.:^(source::QuantumAffine, exponent::Integer)::QuantumAffine = QuantumAffine(source.affine ^ exponent)
Base.:*(quantumaffine::QuantumAffine, region::Subset{Position}) = quantumaffine.affine * region
Base.:*(quantumaffine::QuantumAffine, zone::Subset{Momentum}) = quantumaffine.affine * zone
Base.:*(quantumaffine::QuantumAffine, point::Point) = quantumaffine.affine * point
Spaces.dimension(quantumaffine::QuantumAffine)::Integer = quantumaffine.affine |> dimension

c6 = AffineTransform([cos(π/3) -sin(π/3); sin(π/3) cos(π/3)])

C6 = QuantumAffine(c6)
C6.eigenfunctions

f = BasisFunction(:x => 1, :y => 1im; dimension=2)

functionclass(c6, ct, f)

computecharacterrow(classfunction::BasisFunction) = [(element * classfunction) |> relativephase(classfunction) for element in c6 |> pointgroupelements]
ct * (f |> computecharacterrow)

c6 |> pointgrouporder
eigenfunctions = computeeigenfunctions(c6; dimension=2, functionorderrange=1:3)
(c6 * [eigenfunctions...][12][2]) |> relativephase([eigenfunctions...][12][2])
eigenvalues = [p.first for p in eigenfunctions]
denom = findcomplexdenominator(eigenvalues).denominator
Dict(hashablecomplex(eigenvalue, denom) => (eigenvalue, eigenfunction) for (eigenvalue, eigenfunction) in eigenfunctions) |> values |> collect

triangular = RealSpace([sqrt(3)/2 -1/2; 0. 1.]')

spatialsnappingcalibration((triangular |> origin, triangular & [1/3, 2/3], triangular & [1/2, 1/2]))

c6 = pointgrouptransform([cos(π/3) -sin(π/3); sin(π/3) cos(π/3)]) |> recenter(triangular & [1., 0.])
c6 |> rep
c6.eigenfunctions
findeigenfunction(c6; rankrange=1:3, eigenvalue=1.) |> normalize

p = triangular & [0, 0]

points = Subset(t * p for t in repeat([c6], 6) |> cumprod) + Subset(triangular & [1., 0.])

translate = translation([3, 0], localspace=triangular)
points + translate * points
visualize(points + translate * points, visualspace=euclidean(RealSpace, 2))

c6 * p

triangular |> rep
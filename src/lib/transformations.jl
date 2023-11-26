Base.:(==)(a::BasisFunction, b::BasisFunction) = a.rank == b.rank && (a |> dimension) == (b |> dimension) && isapprox(a |> rep, b |> rep)

Base.:hash(basisfunction::BasisFunction)::UInt = (map(hashablecomplex, basisfunction |> rep), basisfunction |> dimension, basisfunction.rank) |> hash
Base.:isequal(a::BasisFunction, b::BasisFunction)::Bool = a == b

swave::BasisFunction = BasisFunction([1], 0, 0)
export swave

swaveminus::BasisFunction = BasisFunction([-1], 0, 0)
export swaveminus

ppluswave::BasisFunction = BasisFunction([0.5+0.5*sqrt(3)im], 0, 0)
export ppluswave

pminuswave::BasisFunction = BasisFunction([0.5-0.5*sqrt(3)im], 0, 0)
export ppluswave

LinearAlgebra.:normalize(basis::BasisFunction)::BasisFunction = BasisFunction(basis.rep |> normalize, basis.dimension, basis.rank)

Base.:convert(::Type{Vector{Complex}}, source::BasisFunction) = source.rep

Zipper.:dimension(basisfunction::BasisFunction)::Integer = basisfunction.dimension

NAMEINDEXTABLE::Dict{String, Integer} = Dict("x" => 1, "y" => 2, "z" => 3) # The mappings for axis (x y z) in basis functions to tensor indicies.

resolveentrycoords(expression::Symbol)::Vector{Integer} = map(v -> NAMEINDEXTABLE[v], split(expression |> string, ""))

tmatfullindexmappings(dimension::Integer, rank::Integer)::Base.Generator = (
    v |> collect for v in map(reverse, product(repeat([1:dimension], rank)...)))

function tmatcontractrules(dimension::Integer, rank::Integer)::Tuple{OrderedSet{Vector}, Dict{Vector, Set}}
    fullindices = tmatfullindexmappings(dimension, rank)
    ignoreset::Set = Set()
    contracttable::Dict{Vector, Set} = Dict()
    for entry in fullindices
        if entry in ignoreset continue end
        contract::Set = entry |> permutations |> Set
        contracttable[entry] = contract
        union!(ignoreset, Set(v for v in contract if v != entry))
    end
    return OrderedSet(v for v in fullindices if !(v in ignoreset)), contracttable
end

function tmatcontractmatrices(; rank::Integer, dimension::Integer)::Tuple{Matrix, Matrix}
    fulltable::Dict{Vector, Integer} = Dict(v => n for (n, v) in tmatfullindexmappings(dimension, rank) |> enumerate)
    indices::OrderedSet{Vector}, contracttable::Dict{Vector, Set} = tmatcontractrules(dimension, rank)
    invtable::Dict{Integer, Vector} = Dict(n => v for (n, v) in indices |> enumerate)
    contractmatrix::Matrix{Float64} = zeros(Float64, indices |> length, fulltable |> length)
    selectmatrix::Matrix{Float64} = zeros(Float64, indices |> length, fulltable |> length)
    for n in axes(contractmatrix, 1) 
        index::Vector = invtable[n]
        selectmatrix[n, fulltable[index]] = 1
        for contractindex in contracttable[index]
            contractmatrix[n, fulltable[contractindex]] = 1
        end
    end
    return contractmatrix, selectmatrix
end

function Base.:show(io::IO, basisfunction::BasisFunction)
    if basisfunction.rank == 0 # Handles the case of rank 0 basis function.
        print(io, string("$(typeof(basisfunction))(rank=$(basisfunction.rank), $(basisfunction |> rep |> first))"))
        return
    end

    indices::OrderedSet{Vector}, _ = tmatcontractrules(basisfunction.dimension, basisfunction.rank)
    invindexmap::Dict{Integer, Vector} = Dict(n => v for (n, v) in indices |> enumerate)
    indexnametable::Dict{Integer, String} = Dict(i => s for (s, i) in NAMEINDEXTABLE)
    function generatesymbol(info::Tuple)::String
        if isapprox(info |> last |> abs, 0, atol=1e-10) return "" end
        coords::Vector = invindexmap[info |> first]
        printnumber::Number = round(info |> last, digits=4)
        return "($printnumber)*" * reduce(*, Iterators.map(c -> indexnametable[c], coords))
    end
    expression::String = join(filter(v -> v != "", map(generatesymbol, basisfunction |> rep |> enumerate)), " + ")
    print(io, string("$(typeof(basisfunction))(rank=$(basisfunction.rank), $(expression))"))
end

function BasisFunction(expressions::Pair{Symbol, <:Number}...; dimension::Integer)::BasisFunction
    ranks::Tuple = map(expression -> expression |> first |> string |> length, expressions)
    maxrank::Integer = max(ranks...)
    if min(ranks...) != maxrank
        error("Function with mixed order elements is not a valid basis function!")
    end

    indexmap::Dict{Vector, Integer} = Dict(v => n for (n, v) in tmatfullindexmappings(dimension, maxrank) |> enumerate)
    data::Vector{Complex} = zeros(ComplexF64, indexmap |> length)

    for (symbol, value) in expressions
        coords::Vector = symbol |> resolveentrycoords
        data[indexmap[coords]] = value
    end

    contractmatrix, _ = tmatcontractmatrices(dimension=dimension, rank=maxrank)
    return BasisFunction(contractmatrix * data, dimension, maxrank)
end

function Base.:+(a::BasisFunction, b::BasisFunction)::BasisFunction
    @assert(a.dimension == b.dimension)
    @assert(a.rank == b.rank)
    return BasisFunction((a |> rep) + (b |> rep), a.dimension, a.rank)
end

Base.:iszero(basis::BasisFunction)::Bool = basis |> normalize |> rep |> iszero

fullpointgrouprepresentation(irrep::Matrix; rank::Integer = 1)::Matrix = reduce(kron, repeat([irrep], rank))
export fullpointgrouprepresentation

function pointgrouprepresentation(irrep::Matrix; rank::Integer = 1)::Matrix
    contractmatrix, selectormatrix = tmatcontractmatrices(rank=rank, dimension=size(irrep, 1))
    return contractmatrix * fullpointgrouprepresentation(irrep, rank=rank) * selectormatrix'
end

function computeeigenfunctions(pointgroupmatrix::Matrix; functionorderrange::UnitRange = 1:3)::Base.Iterators.Flatten
    dimension::Integer = pointgroupmatrix |> size |> first
    function computeeigenfunctionsatorder(functionorder::Integer)
        matrixatorder::Matrix = pointgrouprepresentation(pointgroupmatrix; rank=functionorder)
        eigenvalues, eigenvectors = matrixatorder |> eigen
        basisfunctions = (eigenvalue => BasisFunction(eigenvectors[:, n], dimension, functionorder) |> normalize for (n, eigenvalue) in eigenvalues |> enumerate)
        return Iterators.filter(p -> !(p.second |> iszero), basisfunctions)
    end

    return (p for functionorder in functionorderrange |> reverse for p in functionorder |> computeeigenfunctionsatorder)
end

function AffineTransform(
    transformmatrix::Matrix, shiftvector::Vector = zeros(Float64, transformmatrix |> size |> first);
    antiunitary::Bool = false,
    localspace::AffineSpace = euclidean(RealSpace, transformmatrix |> size |> first))::AffineTransform

    eigenfunctions = transformmatrix |> computeeigenfunctions
    eigenvaluehashdenominator::Integer = findcomplexdenominator(v for (v, _) in eigenfunctions; denominatorrange=64:128).denominator
    eigenfunctiontable::Dict{Tuple, BasisFunction} = Dict(hashablecomplex(v |> Complex, eigenvaluehashdenominator) => f for (v, f) in eigenfunctions)
    eigenfunctiontable[hashablecomplex(1 + 0im, eigenvaluehashdenominator)] = swave
    eigenfunctiontable[hashablecomplex(-1 + 0im, eigenvaluehashdenominator)] = swaveminus
    eigenfunctiontable[hashablecomplex(0.5+0.5sqrt(3)im, eigenvaluehashdenominator)] = ppluswave
    eigenfunctiontable[hashablecomplex(0.5-0.5sqrt(3)im, eigenvaluehashdenominator)] = pminuswave
    return AffineTransform(localspace, shiftvector, transformmatrix, eigenfunctiontable, eigenvaluehashdenominator, antiunitary)
end

function addeigenfunction!(transform::AffineTransform; eigenvalue::Number, eigenfunction::BasisFunction)
    eigenvaluekey::Tuple = hashablecomplex(eigenvalue |> Complex, transform.eigenvaluehashdenominator)
    trialphase::Complex = relativephase(transform * eigenfunction, eigenfunction)
    hashablecomplex(trialphase |> Complex, transform.eigenvaluehashdenominator) == eigenvaluekey || error("The eigenvalue is not associated with the eigenfunction!")
    transform.eigenfunctions[eigenvaluekey] = eigenfunction
end
export addeigenfunction!

pointgrouptransform(
    pointgroupmatrix::Matrix;
    dimension::Integer = pointgroupmatrix |> size |> first,
    localspace::RealSpace = euclidean(RealSpace, dimension),
    referencepoint::Offset = localspace |> getorigin,
    antiunitary::Bool = false)::AffineTransform = AffineTransform(pointgroupmatrix, transformationshift(pointgroupmatrix, localspace, referencepoint);
    localspace=localspace, antiunitary=antiunitary)
export pointgrouptransform

recenter(transformation::AffineTransform, center::Offset)::AffineTransform = AffineTransform(
    transformation.transformmatrix, transformationshift(transformation.transformmatrix, transformation.localspace, center);
    localspace=transformation |> getspace, antiunitary=transformation.antiunitary)
export recenter

recenter(center::Offset) = transformation -> recenter(transformation, center)

translation(
    shiftvector::Vector;
    dimension::Integer = shiftvector |> length,
    localspace::RealSpace = euclidean(RealSpace, dimension),
    antiunitary::Bool = false)::AffineTransform = AffineTransform(
    Matrix{Float64}(I, dimension, dimension), shiftvector; localspace=localspace, antiunitary=antiunitary)
export translation

identitytransform(dimension::Integer) = AffineTransform(Matrix{Float64}(I, dimension, dimension))
export identitytransform

function affinematrix(transformation::AffineTransform)::Matrix{Float64}
    shiftrow::Vector = vcat(transformation.shiftvector, [1])
    leftcolumns::Matrix = vcat(
        transformation.transformmatrix, zeros(Float64, transformation.transformmatrix |> size |> first, 1) |> transpose)
    return hcat(leftcolumns, shiftrow)
end
export affinematrix

function transformationshift(transformmatrix::Matrix, localspace::AffineSpace, reference::Point)::Vector
    localreference::Point = lineartransform(localspace, reference)
    return (localreference |> vec) - (transformmatrix * (localreference |> vec))
end

Zipper.:getspace(transformation::AffineTransform)::AffineSpace = transformation.localspace

Base.:convert(::Type{Matrix{Float64}}, source::AffineTransform) = source |> affinematrix

function Base.:(==)(a::AffineTransform, b::AffineTransform)::Bool
    localb::AffineTransform = (a |> getspace) * b
    return isapprox(a |> rep, localb |> rep) && a.antiunitary == localb.antiunitary
end

function Base.:*(space::RealSpace, transformation::AffineTransform)::AffineTransform
    if space |> dimension != transformation |> dimension
        error("Dimension mismatch!")
    end
    relativebasis::Matrix = (space |> getbasis |> inv) * (transformation |> getspace |> getbasis)
    transformmatrix::Matrix = relativebasis * (transformation.transformmatrix) * (relativebasis |> inv)
    shiftvector::Vector = lineartransform(space, transformation.shiftvector ∈ (transformation |> getspace)) |> vec
    return AffineTransform(
        transformmatrix, shiftvector;
        localspace=space, antiunitary=transformation.antiunitary)
end

function Base.:*(a::AffineTransform, b::AffineTransform)::AffineTransform
    dim::Integer = a |> dimension
    if dim != b |> dimension
        error("Dimension mismatch!")
    end
    localb::AffineTransform = (a |> getspace) * b
    affinematrix::Matrix = (a |> rep) * (localb |> rep)
    transformmatrix::Matrix = affinematrix[1:dim, 1:dim]
    shiftvector::Vector = affinematrix[1:dim, end]
    antiunitary::Bool = a.antiunitary ⊻ localb.antiunitary
    return AffineTransform(
        transformmatrix, shiftvector;
        localspace=a |> getspace, antiunitary=antiunitary)
end

function Base.:inv(transform::AffineTransform)::AffineTransform
    dim::Integer = transform |> dimension
    affinematrix::Matrix = transform |> rep |> inv
    transformmatrix::Matrix = affinematrix[1:dim, 1:dim]
    shiftvector::Vector = affinematrix[1:dim, end]
    return AffineTransform(
        transformmatrix, shiftvector;
        localspace=transform |> getspace, antiunitary=transform.antiunitary)
end

Base.:*(transformation::AffineTransform, space::RealSpace)::RealSpace = RealSpace((transformation.transformmatrix) * (space |> rep))

function Base.:^(source::AffineTransform, exponent::Number)::AffineTransform
    if exponent == 0
        return AffineTransform(Matrix{Float64}(I, source |> dimension, source |> dimension))
    end
    return reduce(*, repeat([source], exponent))
end

function Base.:*(transformation::AffineTransform, region::Subset{Offset})::Subset{Offset}
    nativetransformation::AffineTransform = (region |> getspace) * transformation
    transformed::Base.Generator = (
        Point(((nativetransformation |> rep) * vcat(point |> vec, [1]))[1:end - 1], point |> getspace)
        for point in region)
    return Subset(transformed)
end

function Base.:*(transformation::AffineTransform, zone::Subset{Momentum})::Subset{Momentum}
    kspace::MomentumSpace = zone |> getspace
    realspace::RealSpace = convert(RealSpace, kspace)
    nativetransformation::AffineTransform = realspace * transformation
    kspacerep::Matrix = Matrix(nativetransformation.transformmatrix |> transpose |> inv)
    return Subset(Point(kspacerep * (k |> vec), kspace) for k in zone)
end

Base.:*(transformation::AffineTransform, point::Point) = (transformation * Subset(point)) |> first

Zipper.:dimension(transformation::AffineTransform)::Integer = transformation.transformmatrix |> size |> first

pointgrouprepresentation(transformation::AffineTransform; rank::Integer = 1)::Matrix = pointgrouprepresentation(transformation.transformmatrix; rank=rank)

function findeigenfunction(transformation::AffineTransform; eigenvalue::Number = 1)::BasisFunction
    eigenvaluekey::Tuple = hashablecomplex(eigenvalue |> Complex, transformation.eigenvaluehashdenominator)
    if !haskey(transformation.eigenfunctions, eigenvaluekey)
        error("No eigenfunction is found for eigenvalue $(eigenvalue)!")
    end
    return transformation.eigenfunctions[eigenvaluekey]
end
export findeigenfunction

function Base.:*(transformation::AffineTransform, basisfunction::BasisFunction)::BasisFunction
    if basisfunction.rank == 0
        return basisfunction
    end

    if dimension(transformation) != dimension(basisfunction)
        error("Dimension mismatch!")
    end

    matrix::Matrix = pointgrouprepresentation(transformation; rank=basisfunction.rank)
    # TODO: Is this really the case?
    functionrep::Vector = transformation.antiunitary ? basisfunction |> rep |> conj : basisfunction |> rep  # Handling anti-unitary.
    transformed::Vector = matrix * functionrep

    return BasisFunction(transformed, basisfunction |> dimension, basisfunction.rank)
end

function pointgroupelements(pointgroup::AffineTransform; maxelements=128)::Vector{AffineTransform}
    identity::AffineTransform = pointgroup ^ 0
    elements::Vector{AffineTransform} = [identity]
    for n in 1:(maxelements - 1)
        current = pointgroup ^ n
        if current == identity
            break
        end
        push!(elements, current)
    end
    return elements
end
export pointgroupelements

pointgrouporder(pointgroup::AffineTransform; maxorder=128)::Integer = pointgroupelements(pointgroup; maxelements=maxorder) |> length
export pointgrouporder

function relativephase(target::BasisFunction, ref::BasisFunction)::Complex
    normalizedref::Vector = ref |> normalize |> rep
    normalizedtarget::Vector = target |> normalize |> rep
    trialphase::Number = dot(normalizedref, normalizedtarget)
    if !isapprox(dot(normalizedtarget ./ trialphase, normalizedref), 1)
        error("Target function is not a gauged version of the reference function.")
    end
    return trialphase
end
export relativephase

relativephase(ref::BasisFunction) = target::BasisFunction -> relativephase(target, ref)

Zipper.:dimension(scale::Scale)::Integer = scale |> rep |> size |> first
Zipper.:getspace(scale::Scale)::RealSpace = scale.localspace

Base.:convert(::Type{Matrix{Float64}}, source::Scale) = source.rep

Base.:inv(scale::Scale)::Scale = Scale(scale |> rep |> inv, scale |> getspace)

Base.:*(a::Scale, b::Scale)::Scale = Scale((a |> rep) * ((a |> getspace) * b |> rep), a |> getspace)

function Base.:*(space::RealSpace, scale::Scale)::Scale
    space |> dimension == scale |> dimension || error("Dimension mismatch!")
    relativebasis::Matrix = (space |> getbasis |> inv) * (scale |> getspace |> getbasis)
    scalematrix::Matrix = relativebasis * (scale |> rep) * (relativebasis |> inv)
    return Scale(scalematrix, space)
end

Base.:*(scale::Scale, space::RealSpace)::RealSpace = RealSpace(*(space * scale |> rep, space |> rep))

Base.:*(scale::Scale, space::MomentumSpace)::MomentumSpace = MomentumSpace(*(convert(RealSpace, space) * scale |> inv |> rep, space |> rep))
Base.:*(scale::Scale, point::Point)::Point = scale * (point |> getspace) * point
Base.:*(scale::Scale, subset::Subset)::Subset = Subset(scale * element for element in subset)

function Base.:*(scale::Scale, crystal::Crystal)::Crystal
    realspace::RealSpace = crystal |> getspace
    snfinput::Matrix{Integer} = map(Integer, vcat(scale |> rep, crystal |> size |> diagm))
    U, S, Vd = snfinput |> dosnf
    snfinput = U[end - dimension(realspace) + 1:end, 1:dimension(realspace)]
    _, bS, bVd = snfinput |> dosnf
    newsizes = bS |> diag
    newrelativebasis = bVd * (S |> diag |> diagm) * Vd |> transpose

    subunitcell = Subset((Vd * (r |> collect)) ∈ realspace for r in Iterators.product((0:(s-1) for s in S |> diag)...))
    newscale = Scale(newrelativebasis, realspace)
    newunitcell = Subset(newscale * (a + b) |> basispoint for (a, b) in Iterators.product(subunitcell, crystal |> getunitcell))
    return Crystal(newunitcell, newsizes)
end

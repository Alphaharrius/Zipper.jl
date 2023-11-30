using Revise, LinearAlgebra, AbstractAlgebra, SmithNormalForm
using Zipper

square = euclidean(RealSpace, 2)
point = [1/2, 1/2] ∈ square
spatialsnappingcalibration([point])

kspace = convert(MomentumSpace, square)

c4 = pointgrouptransform([0 -1; 1 0])

unitcell = Subset(point)
crystal = Crystal(unitcell, [24, 24])
reciprocalhashcalibration(crystal.sizes)

m = quantize(:pos, unitcell, 1) |> first

tₙ = ComplexF64(-1.)
bonds::FockMap = bondmap([
    (m, m |> setattr(:offset => [1, 0] ∈ square)) => tₙ,
    (m, m |> setattr(:offset => [0, 1] ∈ square)) => tₙ])

energyspectrum = computeenergyspectrum(bonds, crystal=crystal)
energyspectrum |> visualize

groundstates::CrystalSpectrum = groundstatespectrum(energyspectrum, perunitcellfillings=0.5)
groundstates |> visualize
groundstateprojector = groundstates |> crystalprojector

C = idmap(groundstateprojector.outspace) - groundstateprojector

correlations = C

correlations |> crystalspectrum |> visualize

crystalfock = correlations |> getoutspace

function dosnf(matrix::Matrix)
    snf = smith(matrix)
    return snf.S, diagm(snf), snf.T
end

    realspace::RealSpace = crystal |> getspace
    snfinput = map(Integer, vcat(scale |> rep, crystal |> size |> diagm))
    U, S, Vd = snfinput |> computesnf
    snfinput = U[end - dimension(realspace) + 1:end, 1:dimension(realspace)]
    bU, bS, bVd = snfinput |> computesnf
    newsizes = bS |> diag
    newrelativebasis = bVd * (S |> diag |> diagm) * Vd |> transpose

    subunitcell = Subset((Vd * (r |> collect)) ∈ realspace for r in Iterators.product((0:(s-1) for s in S |> diag)...))
    newscale = Scale(newrelativebasis, realspace)
    newunitcell = Subset(newscale * (a + b) |> basispoint for (a, b) in Iterators.product(subunitcell, crystal |> getunitcell))
    newcrystal = Crystal(newunitcell, newsizes)

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

function Base.:*(scale::Scale, crystal::Crystal)::Crystal
    realspace::RealSpace = crystal |> getspace
    snfinput::Matrix{Integer} = map(Integer, vcat(scale |> rep, crystal |> size |> diagm))
    U, S, Vd = snfinput |> computesnf
    snfinput = U[end - dimension(realspace) + 1:end, 1:dimension(realspace)]
    _, bS, bVd = snfinput |> computesnf
    newsizes = bS |> diag
    newrelativebasis = bVd * (S |> diag |> diagm) * Vd |> transpose

    subunitcell = Subset((Vd * (r |> collect)) ∈ realspace for r in Iterators.product((0:(s-1) for s in S |> diag)...))
    newscale = Scale(newrelativebasis, realspace)
    newunitcell = Subset(newscale * (a + b) for (a, b) in Iterators.product(subunitcell, crystal |> getunitcell))
    return Crystal(newunitcell, newsizes)
end

scale = Scale([2 0; 0 2], square)
scale * crystal |> getspace |> rep
blockresult = blocking(:action => scale, :correlations => correlations, :crystal => crystalfock |> getcrystal)

H = energyspectrum |> FockMap
blockresult[:transformer] * H * blockresult[:transformer]' |> crystalspectrum |> visualize

circleregion = getsphericalregion(from=square |> getorigin, generators=Subset(square |> getbasisvectors), symmetries=[c4], radius=10, metricspace=square)
circleregion |> visualize

normaldirection = [1, 1] ∈ square
orthodirection = c4 * normaldirection

crosssectionfilter(point::Point) = 0 < dot(point, normaldirection |> normalize) < norm(normaldirection) && (dot(point, orthodirection) |> abs) < 4

blockedcrystal = blockresult[:crystal]
blockedcrystal |> getspace |> rep

blockedcrystal |> sitepoints |> visualize

function extendedcrystalrestriction(crystal::Crystal, generatingvector::Offset)
    snfinput::Matrix{Integer} = map(Integer, vcat(generatingvector |> vec |> transpose, crystal |> size |> diagm)) 
    _, S, Vd = snfinput |> dosnf
    scale::Scale = Scale((S |> diag |> diagm) * Vd, crystal |> getspace)
    show(scale |> rep)
    barecrystal::Crystal = Crystal(crystal |> getspace |> getorigin |> Subset, crystal |> size)
    return scale * barecrystal
end

embeddedcrystal = extendedcrystalrestriction(blockedcrystal, [1, 1] ∈ square)
embeddedcrystal |> getspace |> getbasisvectors |> collect

embeddedcrystal |> sitepoints |> visualize


function Base.:*(scale::Scale, crystal::Crystal)::Crystal
    oldspace::RealSpace = crystal |> getspace
    snf = vcat(scale |> rep, crystal.sizes |> diagm) |> smith
    boundarysnf = snf.S[end - dimension(oldspace) + 1:end, 1:dimension(oldspace)] |> smith
    Δ::Matrix{Float64} = snf |> diagm
    newbasiscoords::Matrix{Float64} = boundarysnf.T * (Δ |> diag |> diagm) * snf.T
    blockingpoints::Base.Generator = (Point(collect(coord)' * snf.T, oldspace) for coord in Iterators.product((0:size - 1 for size in diag(Δ))...))
    relativescale::Scale = Scale(newbasiscoords, crystal |> getspace)
    scaledunitcell::Subset{Offset} = Subset(relativescale * (a + b) for (a, b) in Iterators.product(blockingpoints, crystal.unitcell))
    return Crystal(scaledunitcell, diag(diagm(boundarysnf)))
end

    oldspace::RealSpace = crystal |> getspace
    snf = vcat(scale |> rep, crystal.sizes |> diagm) |> smith
    U, S, Vd = map(Integer, vcat(scale |> rep, crystal.sizes |> diagm)) |> computesnf
    boundarysnf = snf.S[end - dimension(oldspace) + 1:end, 1:dimension(oldspace)] |> smith
    Δ::Matrix{Float64} = snf |> diagm
    newbasiscoords::Matrix{Float64} = boundarysnf.T * (Δ |> diag |> diagm) * snf.T
    blockingpoints::Base.Generator = (Point(collect(coord)' * snf.T, oldspace) for coord in Iterators.product((0:size - 1 for size in diag(Δ))...))
    relativescale::Scale = Scale(newbasiscoords, crystal |> getspace)
    scaledunitcell::Subset{Offset} = Subset(relativescale * (a + b) for (a, b) in Iterators.product(blockingpoints, crystal.unitcell))
    return Crystal(scaledunitcell, diag(diagm(boundarysnf)))

function Base.:*(scale::Scale, crystal::Crystal)::Crystal
    realspace::RealSpace = crystal |> getspace
    snfinput::Matrix{Integer} = map(Integer, vcat(scale |> rep, crystal |> size |> diagm))
    U, S, Vd = snfinput |> computesnf
    snfinput = U[end - dimension(realspace) + 1:end, 1:dimension(realspace)]
    _, bS, bVd = snfinput |> computesnf
    newsizes = bS |> diag
    newrelativebasis = bVd * (S |> diag |> diagm) * Vd
    show(newrelativebasis)

    subunitcell = Subset((Vd * (r |> collect)) ∈ realspace for r in Iterators.product((0:(s-1) for s in S |> diag)...))
    newscale = Scale(newrelativebasis, realspace)
    newunitcell = Subset(newscale * (a + b) for (a, b) in Iterators.product(subunitcell, crystal |> getunitcell))
    return Crystal(newunitcell, newsizes)
end

blockedcrystal = blockresult[:crystal]

blockedcrystal |> getspace |> rep
normaldirection
normaldirection = [1, 1] ∈ square

function extendedcrystalrestriction(crystal::Crystal, generatingvector::Offset)
    snfinput::Matrix{Integer} = map(Integer, vcat(generatingvector |> vec |> transpose, crystal |> size |> diagm)) 
    _, S, Vd = snfinput |> computesnf
    scale::Scale = Scale((S |> diag |> diagm) * Vd, crystal |> getspace)
    barecrystal::Crystal = Crystal(crystal |> getspace |> getorigin |> Subset, crystal |> size)

    return scale * barecrystal
end

snfinput = map(Integer, vcat(normaldirection |> vec |> transpose, blockedcrystal |> size |> diagm))
U, S, Vd = snfinput |> dosnf
U * S * Vd
scalematrix = (S |> diag |> diagm) * Vd
scale = Scale(scalematrix, blockedcrystal |> getspace)
barecrystal = Crystal(blockedcrystal |> getspace |> getorigin |> Subset, blockedcrystal |> size)
barecrystal |> getspace |> rep

scale |> rep

function Base.:*(scale::Scale, crystal::Crystal)::Crystal
    realspace::RealSpace = crystal |> getspace
    snfinput::Matrix{Integer} = map(Integer, vcat(scale |> rep, crystal |> size |> diagm))
    U, S, Vd = snfinput |> computesnf
    snfinput = U[end - dimension(realspace) + 1:end, 1:dimension(realspace)]
    _, bS, bVd = snfinput |> computesnf
    newsizes = bS |> diag
    newrelativebasis = bVd * (S |> diag |> diagm) * Vd
    show(newrelativebasis)

    subunitcell = Subset((Vd * (r |> collect)) ∈ realspace for r in Iterators.product((0:(s-1) for s in S |> diag)...))
    newscale = Scale(newrelativebasis, realspace)
    newunitcell = Subset(newscale * (a + b) for (a, b) in Iterators.product(subunitcell, crystal |> getunitcell))
    return Crystal(newunitcell, newsizes)
end

newcrystal = scale * barecrystal

space = barecrystal |> getspace
snfinput = map(Integer, vcat(scale |> rep, barecrystal |> size |> diagm))
U, S, Vd = snfinput |> computesnf
snfinput = U[end - dimension(space) + 1:end, 1:dimension(space)]
bU, bS, bVd = snfinput |> computesnf
newsizes = bS |> diag
newrelativebasis = bVd * (S |> diag |> diagm) * Vd |> transpose

subunitcell = Subset((Vd * (r |> collect)) ∈ space for r in Iterators.product((0:(s-1) for s in S |> diag)...))
newscale = Scale(newrelativebasis, space)
newunitcell = Subset(newscale * (a + b) for (a, b) in Iterators.product(subunitcell, barecrystal |> getunitcell))
newunitcell |> collect
newcrystal = Crystal(newunitcell, newsizes)





visualize(newcrystal |> sitepoints, barecrystal |> sitepoints)

visualize(circleregion, filter(crosssectionfilter, circleregion))
map(Integer, vcat(normaldirection |> vec |> transpose, blockedcrystal |> size |> diagm))
S, T, U = matrix(ZZ, map(Integer, vcat(normaldirection |> vec |> transpose, blockedcrystal |> size |> diagm))) |> snf_with_transform
sU = inv(T)
sVd = inv(U)
S

sU * S * sVd

(S |> Matrix |> diag |> diagm) * Matrix(U |> inv)

snf = vcat(normaldirection |> vec |> transpose, blockedcrystal |> size |> diagm) |> smith
extendscale = Scale((snf |> diagm |> diag |> diagm) * snf.T, blockedcrystal |> getspace)
barecrystal = Crystal(blockedcrystal |> getspace |> getorigin |> Subset, blockedcrystal |> size)

extendscale |> rep
space::RealSpace = barecrystal |> getspace

S, invU, invVd = matrix(ZZ, map(Integer, vcat(extendscale |> rep, barecrystal |> size |> diagm))) |> snf_with_transform
newsize, invboundaryU, invboundaryVd = inv(invU)[end - dimension(space) + 1:end, 1:dimension(space)] |> snf_with_transform
newbasis = map(Int64, inv(invboundaryVd) * matrix(ZZ, S |> Matrix |> diag |> diagm) * inv(invVd) |> Matrix)
Vd = map(Int64, invVd |> inv |> Matrix)
newpoints = (Point((shift |> collect |> transpose) * Vd, space) for shift in product((0:size - 1 for size in map(Int64, S |> Matrix |> diag))...))
newscale = Scale(reverse(newbasis, dims=2), space)
newunitcell = Subset(newscale * (a + b) for (a, b) in Iterators.product(newpoints, barecrystal |> getunitcell))
newcrystal = Crystal(newunitcell, newsize |> Matrix |> diag)

newcrystal |> getunitcell |> collect
newcrystal |> sitepoints |> visualize


barecrystal = Crystal(blockedcrystal |> getspace |> getorigin |> Subset, blockedcrystal |> size)
extendscale |> rep
scaledcrystal = extendscale * barecrystal
scaledcrystal |> getspace |> getbasisvectors |> collect
for (basis, d) in zip(scaledcrystal |> getspace |> rep, scaledcrystal |> size)
end
scaledcrystal |> getspace |> rep
scaledcrystal |> getspace |> getbasisvectors |> collect
fr = Iterators.filter(zip(scaledcrystal |> getspace |> getbasisvectors, scaledcrystal |> size)) do v
    return v |> last > 1
end
fr |> collect |> first |> first |> getspace |> rep
affinespace(fr |> collect |> first, )

visualize(barecrystal |> sitepoints, scaledcrystal |> sitepoints)

scaledcrystal |> getspace |> rep
extendscale |> rep
visualize(barecrystal |> sitepoints, scaledcrystal |> sitepoints)

extendunitcell = extendscale * blockedcrystal |> getunitcell
extendcrystal = Crystal(Subset(extendunitcell[2], extendunitcell[1]), [3, 1])
extendcrystal |> getspace |> getbasisvectors |> collect

extendcrystal |> sitepoints |> visualize

snf = vcat(extendscale |> rep, blockedcrystal |> size |> diagm) |> smith
snf.S
snf |> diagm |> diag
snf.T

barecrystal = Crystal(blockedcrystal |> getspace |> getorigin |> Subset, blockedcrystal |> size)
barecrystal |> sitepoints |> visualize
extendscale * barecrystal |> getunitcell |> visualize


bS, bT, bU = inv(snf.T)[end - dimension(blockedcrystal |> getspace) + 1:end, 1:dimension(blockedcrystal |> getspace)] |> snf_with_transform
Matrix(U |> inv) * (S |> Matrix |> diag |> diagm) * Matrix(bU)

bsnf.S
A = bsnf.T * (snf |> diagm |> diag |> diagm) * snf.T
Scale(A, blockedcrystal |> getspace) * point == point


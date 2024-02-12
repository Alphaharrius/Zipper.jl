using Base.Iterators
using Zipper

square = euclidean(RealSpace, 2)
point = [1/2, 1/2] ∈ square
spatialsnappingcalibration([point])

kspace = convert(MomentumSpace, square)

c4 = pointgrouptransform([0 -1; 1 0])

unitcell = Subset(point)
crystal = Crystal(unitcell, [384, 384])
reciprocalhashcalibration(crystal.sizes)

using OrderedCollections
hash([Subset([1/2, 1/2] ∈ square), [384, 384]])

bz = brillouinzone(crystal)

getindexk(crystal::Crystal, k::Momentum) = vec(k) .* size(crystal)

kindexoperator(crystal::Crystal) = [1, (size(crystal)[1:end-1]|>cumprod)...]

Base.getindex(crystal::Crystal, k::Momentum)::Integer = kindexoperator(crystal)' * getindexk(crystal, k) + 1

Base.:in(k::Momentum, crystal::Crystal) = getspace(k) == convert(MomentumSpace, getspace(crystal)) && getindexk(crystal, k)|>sum|>isinteger

k = bz[129]
k ∈ crystal
kindex(crystal, bz[123214])

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

crystalpoints::Subset{Point} = latticepoints(newcrystal)
newmodes::Subset{Mode} = removeattr(rep(blockedcorrelation.outspace)[1], :offset)
physicalmodes::Subset{Mode} = spanoffset(newmodes, crystalpoints)
restrictedregion::Subset{Mode} = filter(circularfilter(origin(euclidean(RealSpace, 2)), 1.0), physicalmodes)
restrictedfock::FockSpace = FockSpace(restrictedregion)

restrictfourier::FockMap = fourier(blockedcorrelation.inspace, restrictedfock) / sqrt(vol(newcrystal))
visualize(restrictfourier, rowrange=1:32)
restrictedcorrelation::FockMap = restrictfourier' * blockedcorrelation * restrictfourier
visualize(restrictedcorrelation)
crvs = eigvalsh(restrictedcorrelation)
plot(scatter(y=map(p -> p.second, crvs), mode="markers"))

restrictedunitary::FockMap = eigvecsh(restrictedcorrelation)
visualize(restrictedunitary)

emode1::Mode = orderedmodes(restrictedunitary.inspace)[4]
emode2::Mode = orderedmodes(restrictedunitary.inspace)[5]
mr1::FockMap = columns(restrictedunitary, FockSpace(Subset([emode1])))
mr2::FockMap = columns(restrictedunitary, FockSpace(Subset([emode2])))
mr = FockMap(mr1.outspace, mr1.inspace, rep(mr1) + 1im * rep(mr2))
values = columnspec(mr)

visualize(values)

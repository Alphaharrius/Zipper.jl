include("../src/spaces.jl")
include("../src/geometries.jl")
include("../src/transformations.jl")
include("../src/quantum.jl")
include("../src/physical.jl")
include("../src/plotting.jl")
include("../src/quantumtransformations.jl")
include("../src/renormalization.jl")

using PlotlyJS, SmithNormalForm, LinearAlgebra, OrderedCollections, SparseArrays, Combinatorics
using ..Spaces, ..Geometries, ..Quantum, ..Transformations, ..Plotting, ..QuantumTransformations, ..Physical, ..Renormalization

triangular = RealSpace([sqrt(3)/2 -1/2; 0. 1.]')
kspace = convert(MomentumSpace, triangular)

c6 = triangular * isometrictransformation([cos(π/3) -sin(π/3); sin(π/3) cos(π/3)])

unitcell = Subset(triangular & [1/3, 2/3], triangular & [2/3, 1/3])
crystal = Crystal(unitcell, [64, 64])

modes::Subset{Mode} = quantize(:pos, unitcell, 1)
m0, m1 = members(modes)

tₙ = ComplexF64(-1.)
bonds::FockMap = bondmap([
    (m0, m1) => tₙ,
    (m0, setattr(m1, :offset => Point([-1, 0], triangular))) => tₙ,
    (m0, setattr(m1, :offset => Point([0, 1], triangular))) => tₙ])

𝐻::FockMap = hamiltonian(crystal, bonds)
𝐶::FockMap = groundstatecorrelations(𝐻)

crystalspec = crystalspectrum(𝐻)

corrspec = 𝐶 |> crystalspectrum

visualize(crystalspec, title="Physical Hamiltonian")
visualize(corrspec, title="Physical Correlation")

crystalfock = 𝐻.outspace

blockresult = blocking(:action => Scale([2 0; 0 2]), :correlations => 𝐶, :crystal => crystal)

newcrystal = blockresult[:crystal]

blockedH = blockresult[:transformer] * 𝐻 * blockresult[:transformer]'

blockedHspectrum = blockedH |> crystalspectrum
visualize(blockedHspectrum, title="Blocked Hamiltonian")

crystalpoints::Subset{Position} = latticepoints(newcrystal)
newmodes::Subset{Mode} = quantize(:pos, newcrystal.unitcell, 1)
physicalmodes::Subset{Mode} = spanoffset(newmodes, crystalpoints)
restrictedregion::Subset{Mode} = filter(circularfilter(origin(euclidean(RealSpace, 2)), 2.0), physicalmodes)
restrictedfock::FockSpace = FockSpace(restrictedregion)

dH = globaldistillerhamiltonian(;
    correlations=blockresult[:correlations],
    restrictspace=restrictedfock,
    localisometryselectionstrategy=frozenselectionbycount(3),
    symmetry=c6)

gdspectrum = dH |> crystalspectrum

gdspectrum |> visualize

distillresult = dH |> distillation

seed = [regionalwannierseeding(blockresult[:correlations], restrictedfock, symmetry=c6)...][1]

function wannierprojection(; crystalisometries::Dict{Momentum, FockMap}, crystal::Crystal, seed::FockMap, svdorthothreshold::Number = 1e-1)
    wannierunitcell::Subset{Position} = Subset(mode |> getattr(:pos) for mode in seed.inspace |> orderedmodes)
    wanniercrystal::Crystal = Crystal(wannierunitcell, crystal.sizes)
    overlaps = ((k, isometry, isometry' * rows(seed, isometry.outspace)) for (k, isometry) in crystalisometries)
    function approximateisometry(k::Momentum, isometry::FockMap, overlap::FockMap)::FockMap
        U, Σ, Vt = overlap |> svd
        minsvdvalue::Number = minimum(v for (_, v) in Σ)
        if minsvdvalue < svdorthothreshold
            # @warn("Precarious Wannierization with a minimum overlap svd value of $(minsvdvalue)!")
        end
        unitary::FockMap = U * Vt
        approximated::FockMap = isometry * unitary

        return FockMap(approximated, inspace=FockSpace(approximated.inspace |> orderedmodes |> setattr(:offset => k)), performpermute=false)
    end
    wannierisometry::FockMap = directsum(approximateisometry(k, isometry, overlap) for (k, isometry, overlap) in overlaps)
    inspace::FockSpace = FockSpace(wannierisometry.inspace, reflected=wanniercrystal)
    outspace::FockSpace = FockSpace(wannierisometry.outspace, reflected=crystal)
    return FockMap(wannierisometry, inspace=inspace, outspace=outspace, performpermute=false)
end

filledwannierisometry = wannierprojection(crystalisometries=distillresult[:filled].eigenvectors, crystal=newcrystal, seed=seed)

filledH = filledwannierisometry' * blockedH * filledwannierisometry

filledH |> crystalspectrum |> visualize

emptywannierisometry = wannierprojection(crystalisometries=distillresult[:empty].eigenvectors, crystal=newcrystal, seed=seed)
emptyH = emptywannierisometry' * blockedH * emptywannierisometry
emptyH |> crystalspectrum |> visualize

bigregion::Subset{Mode} = filter(circularfilter(origin(euclidean(RealSpace, 2)), 2.0), physicalmodes)
bigregionfock::FockSpace = FockSpace(bigregion)
Ft = fourier(blockedH.outspace, bigregionfock)

filledisometries = Ft' * filledwannierisometry
columns(filledisometries, (filledisometries.inspace |> orderedmodes)[1] |> Subset |> FockSpace) |> columnspec |> visualize

overlaps = [unitary' * rows(seed, unitary.outspace) for (k, unitary) in distillresult[:filled].eigenvectors]
u, s, vt = overlaps[1] |> svd
[s...]

visualize(u * vt)

seedoutspace = FockSpace(seed.outspace, reflected=blockresult[:crystal])
checker = (c6 * seedoutspace) * seed - seed
visualize(checker, rowrange=1:128)

restrictedregion::Subset{Mode} = filter(circularfilter(origin(euclidean(RealSpace, 2)), 2.), physicalmodes)
restrictedfock::FockSpace = FockSpace(restrictedregion)

localcorr = Zer.regioncorrelations(blockresult[:correlations], restrictedfock)
localcorrspec = localcorr |> eigenspectrum
validgroups = [Iterators.filter(p -> p.first <= 1e-2, localcorrspec |> groupbyeigenvalues)...]

validgroups[1]

localiso = columns(localcorrspec.eigenvectors, validgroups[2].second |> FockSpace)
visualize(localiso)



region = Subset(m |> pos for m in restrictedregion)
visualize(region, visualspace=euclidean(RealSpace, 2))
c6 * region == region

c6T = c6 * localiso.outspace
visualize(c6T)

U = localiso' * c6T * localiso
visualize(U)

function LinearAlgebra.log(fockmap::FockMap)::FockMap
    mat::SparseMatrixCSC = fockmap |> rep |> Matrix |> log |> SparseMatrixCSC
    return FockMap(fockmap.outspace, fockmap.inspace, mat)
end

projectedgenerator = (U |> log) * 1im
hermitiangenerator = (projectedgenerator + projectedgenerator') / 2
unitary = eigvecsh(hermitiangenerator)

transformediso = localiso * unitary
visualize(transformediso' * transformediso)

function symmetrizeseed(seedisometry::FockMap)::FockMap
    symmetrizers = (transformation * seedisometry.outspace for transformation in c6 |> pointgroupelements)
    return reduce(+, symmetrizer * seedisometry for symmetrizer in symmetrizers)
end

symmetrizediso = transformediso
columns(symmetrizediso, (symmetrizediso.inspace |> orderedmodes)[2] |> Subset |> FockSpace) |> columnspec |> visualize

columns(localiso, (localiso.inspace |> orderedmodes)[1] |> Subset |> FockSpace) |> columnspec |> visualize

function symmetrizeseed(seedisometry::FockMap)::FockMap
    symmetrizers = (transformation * seedisometry.outspace for transformation in c6 |> pointgroupelements)
    return reduce(+, symmetrizer * seedisometry for symmetrizer in symmetrizers)
end

slocaliso = localiso |> symmetrizeseed
visualize(slocaliso)

visualize(c6 * slocaliso.outspace)

newcrystalfock = blockresult[:correlations].inspace
F = fourier(newcrystalfock, slocaliso.outspace)

g = F * localiso
visualize(g' * g)

globalisometries::Dict{Momentum, FockMap} = crystalisometries(localisometry=slocaliso, crystalfock=blockresult[:correlations].inspace)
globalisometry::FockMap = directsum(v for (_, v) in globalisometries)
visualize(globalisometry' * globalisometry)

globalisos = crystalisometries(localisometry=localiso, crystalfock=blockresult[:correlations].inspace)
globaliso = FockMap(directsum(v for (_, v) in globalisos), outspace=blockresult[:correlations].inspace)
pseudoiden = globaliso' * globaliso
symmetrizers = [(transformation * globaliso.outspace for transformation in c6 |> pointgroupelements)...]
result = reduce(+, symmetrizer * globaliso for symmetrizer in symmetrizers)
visualize(result)

((c6 * localiso.outspace) * localiso) + localiso

newregion::Subset{Mode} = filter(circularfilter(origin(euclidean(RealSpace, 2)), 8.0), physicalmodes)
newregionfock::FockSpace = FockSpace(newregion)

ifft = Quantum.fourier(blockresult[:correlations].inspace, newregionfock)'

unblockedvec = blockresult[:transformer]' * columns(result, [(result.inspace |> orderedmodes)[1]] |> Subset |> FockSpace)
physfock = Quantum.spanoffset(modes, latticepoints(crystal)) |> FockSpace
ifft = Quantum.fourier(crystalfock, physfock)'
vvec = ifft * FockMap(unblockedvec, outspace=ifft.inspace)
visualize(vvec |> columnspec)


vec = ifft * columns(result, [(result.inspace |> orderedmodes)[1]] |> Subset |> FockSpace)
visualize(vec)
visualize(Quantum.columnspec(vec))

struct BlockDiagEigenSpectrum{T}
    groupedeigenmodes::Dict{T, Vector{Mode}}
    eigenvalues::Dict{Mode, Number}
    eigenvectors::Dict{T, FockMap}
end

function blockdiageigenspectrum(fockmap::FockMap, groupkey::Symbol; renamekey::Symbol = :groupkey)::BlockDiagEigenSpectrum
    if !hassamespan(fockmap.inspace, fockmap.outspace) @error("Not a Hermitian!") end
    hermitian::FockMap = Quantum.permute(fockmap, outspace=fockmap.inspace) # Make sure the inspace and outspace have the same order, with the inspace as reference.

    submaps::Base.Generator = (commonattr(subspace |> orderedmodes, groupkey) => restrict(fockmap, subspace, subspace) for subspace in hermitian.inspace |> subspaces)
    T::Type = submaps |> first |> first |> typeof # Get the type of the grouping attribute.
    eigens::Base.Generator = (key => eigh(submap, renamekey => key) for (key, submap) in submaps)
    groupedeigenvalues::Base.Generator = (key => evals for (key, (evals, _)) in eigens)
    eigenvectors::Dict{T, FockMap} = Dict(key => evecs for (key, (_, evecs)) in eigens)
    groupedeigenmodes::Dict{T, Vector{Mode}} = Dict(key => map(v -> v.first, evals) for (key, evals) in groupedeigenvalues)
    eigenvalues::Dict{Mode, Number} = Dict(mode => v for (_, evals) in groupedeigenvalues for (mode, v) in evals)

    return BlockDiagEigenSpectrum{T}(groupedeigenmodes, eigenvalues, eigenvectors)
end

pseudospec = blockdiageigenspectrum(pseudoiden, :offset; renamekey=(:offset))
pseudospec.eigenvalues

seed = [regionalwannierseeding(blockresult[:correlations], restrictedfock, symmetry=c6)...][1]
visualize(seed, rowrange=1:1024)
5e-2

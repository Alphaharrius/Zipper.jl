include("../src/spaces.jl")
include("../src/geometries.jl")
include("../src/transformations.jl")
include("../src/quantum.jl")
include("../src/physical.jl")
include("../src/plotting.jl")
include("../src/quantumtransformations.jl")
include("../src/renormalization.jl")

using Plotly, SmithNormalForm, LinearAlgebra, OrderedCollections, SparseArrays, Combinatorics
using ..Spaces, ..Geometries, ..Quantum, ..Transformations, ..Plotting, ..QuantumTransformations, ..Physical, ..Renormalization

triangular = RealSpace([sqrt(3)/2 -1/2; 0. 1.]')
kspace = convert(MomentumSpace, triangular)

pa = triangular & [1/3, 2/3]
pb = triangular & [2/3, 1/3]
pc = (pa + pb) / 2
spatialsnappingcalibration((pa, pb, pc))

c6 = pointgrouptransform([cos(π/3) -sin(π/3); sin(π/3) cos(π/3)])

unitcell = Subset(pa, pb)
crystal = Crystal(unitcell, [32, 32])
reciprocalhashcalibration(crystal.sizes)

modes::Subset{Mode} = quantize(:pos, unitcell, 1)
m0, m1 = members(modes)

tₙ = ComplexF64(-1.)
bonds::FockMap = bondmap([
    (m0, m1) => tₙ,
    (m0, setattr(m1, :offset => Point([-1, 0], triangular))) => tₙ,
    (m0, setattr(m1, :offset => Point([0, 1], triangular))) => tₙ])

H::FockMap = hamiltonian(crystal, bonds)
C::FockMap = groundstatecorrelations(H)

energyspectrum = crystalspectrum(H)
correlationspectrum = C |> crystalspectrum

visualize(energyspectrum, title="Physical Hamiltonian")
visualize(correlationspectrum, title="Physical Correlation")

crystalfock = H.outspace

scale = Scale([2 0; 0 2])
blockresult = blocking(:action => scale, :correlations => C, :crystal => crystal)

newcrystal = blockresult[:crystal]

blockedH = blockresult[:transformer] * H * blockresult[:transformer]'

blockedHspectrum = blockedH |> crystalspectrum
visualize(blockedHspectrum, title="Blocked Hamiltonian")

commutation(a::FockMap, b::FockMap)::FockMap = (a * b - b * a)

commute(c6 * blockedH.outspace, blockedH; eps=5e-7)

function circularregionmodes(origin::Offset, physicalmodes::Subset{Mode}, radius::Number)::Subset{Mode}
    physicalnorm = m -> m |> pos |> euclidean |> norm
    return filter(p -> physicalnorm(p - origin) < radius, physicalmodes)
end

crystalpoints::Subset{Offset} = latticepoints(newcrystal)
newmodes::Subset{Mode} = quantize(:pos, newcrystal.unitcell, 1)
physicalmodes::Subset{Mode} = spanoffset(newmodes, crystalpoints)

frozenseedingregion::Subset{Mode} = circularregionmodes(triangular |> origin, physicalmodes, 2.0)
frozenseedingfock::FockSpace = FockSpace(frozenseedingregion)
visualize(Subset(m |> pos |> euclidean for m in frozenseedingregion), visualspace=euclidean(RealSpace, 2))

courierseedingcenter::Offset = (newcrystal |> getspace) & [2/3, 1/3]
courierseedingregion::Subset{Mode} = circularregionmodes(courierseedingcenter, physicalmodes, 1.8)
courierseedingfock::FockSpace = FockSpace(courierseedingregion)
visualize(Subset(m |> pos |> euclidean for m in courierseedingregion), visualspace=euclidean(RealSpace, 2))

dH = globaldistillerhamiltonian(;
    correlations=blockresult[:correlations],
    restrictspace=frozenseedingfock,
    localisometryselectionstrategy=frozenselectionbycount(3),
    symmetry=c6)

gdspectrum = dH |> crystalspectrum
gdspectrum |> visualize

distillresult = dH |> distillation

filledprojector = distillresult[:filled] |> crystalprojector
filledcorrelation = idmap(filledprojector.outspace, filledprojector.outspace) - filledprojector
filledcorrelation |> crystalspectrum |> visualize

filledseed = [regionalwannierseeding(filledcorrelation, frozenseedingfock, symmetry=c6, seedsgroupingprecision=1e-3)...][1]

c3 = (c6^2) |> recenter(courierseedingcenter)

using Statistics
function Transformations.:relativephase(a::FockMap, b::FockMap)::Complex
    elementwisediv::FockMap = FockMap(a.outspace, a.inspace, (a |> rep) ./ (b |> Quantum.permute(outspace=a.outspace, inspace=a.inspace) |> rep))
    trialphase::Complex = (elementwisediv |> maximum, elementwisediv |> minimum) |> mean
    trialphase |> println
    (a / trialphase - b) |> makezero |> iszero || error("The parameters aren't equivalent up to a phase!")
    return trialphase
end

courierprojector = distillresult[:courier] |> crystalprojector
couriercorrelation = idmap(courierprojector.outspace, courierprojector.outspace) - courierprojector

commutation(c6 * couriercorrelation.outspace, couriercorrelation) |> maximum



courierseed = [regionalwannierseeding(couriercorrelation, courierseedingfock, symmetry=c3)...][1]
fullcourierseed = courierseed + (c6 * courierseed.outspace) * courierseed * (c6 * courierseed.inspace)'

fullcourierseed.inspace |> modeattrs

c6rep = c6 * courierseed.outspace

symmetrizedseed = reduce(+, (t * fullcourierseed.outspace) * fullcourierseed for t in c6 |> pointgroupelements)

function Renormalization.crystalisometry(; localisometry::FockMap, crystalfock::FockSpace{Crystal})::FockMap
    isometries::Dict{Momentum, FockMap} = crystalisometries(
        localisometry=localisometry, crystalfock=crystalfock, addinspacemomentuminfo=true)
    isometryunitcell::Subset{Offset} = Subset(mode |> getattr(:pos) for mode in localisometry.inspace |> orderedmodes)
    isometrycrystal::Crystal = Crystal(isometryunitcell, (crystalfock |> getcrystal).sizes)
    isometry::FockMap = (isometry for (_, isometry) in isometries) |> directsum
    isometrycrystalfock::CrystalFock = FockSpace(isometry.inspace, reflected=isometrycrystal)
    return FockMap(isometry, inspace=isometrycrystalfock, outspace=crystalfock)
end

crystalseed = crystalisometry(localisometry=fullcourierseed, crystalfock=blockedH.outspace)

(c6 * crystalseed.outspace) * crystalseed * (c6 * crystalseed.inspace)' - crystalseed |> maximum

(c6 * (crystalseed.inspace |> unitcellfock)).outspace |> modeattrs

c3o = c6 ^ 2
crystalcourierseed - (c3o * crystalcourierseed.outspace) * crystalcourierseed * (c3o * crystalcourierseed.inspace)' |> maximum

courierseed + courierseed2

courierseed - (c3 * courierseed.outspace) * courierseed |> maximum
courierseed.inspace |> modeattrs
courierseed2 = (c6 * courierseed.outspace) * courierseed * (c6 * courierseed.inspace)'

courierseedregion = Subset(m |> pos for m in courierseed.outspace |> orderedmodes)
courierseedregion2 = Subset(m |> pos for m in courierseed2.outspace |> orderedmodes)
couriercenter2 = courierseedregion2 |> center
c3_2 = c6^2 |> recenter(couriercenter2)
c3_2 * courierseedregion2 == courierseedregion2
c3_2 * courierseed2.outspace * courierseed2 - courierseed2 |> maximum

fullcourierseed = courierseed + courierseed2
fullcourierregion = Subset(m |> pos for m in fullcourierseed.outspace |> orderedmodes)
visualize(fullcourierregion, c6 * fullcourierregion, visualspace=euclidean(RealSpace, 2))

crystalcourierseeds = crystalisometries(localisometry=fullcourierseed, crystalfock=blockedH.outspace, addinspacemomentuminfo=false)
crystalcourierseed = FockMap((seed for (_, seed) in crystalcourierseeds) |> directsum, outspace=blockedH.outspace)

visualize(crystalcourierseed' * crystalcourierseed)
visualize(filledseed' * filledseed)

(c6 * crystalcourierseed.outspace) * crystalcourierseed - crystalcourierseed |> maximum

commutation(c6 * crystalcourierseed.outspace, crystalcourierseed * crystalcourierseed') |> maximum

visualize(courierseedregion, courierseedregion2, visualspace=euclidean(RealSpace, 2))

courierseedproj = crystalcourierseed * crystalcourierseed'

wanniercourier = wannierprojection(crystalisometries=distillresult[:courier].eigenvectors, crystal=newcrystal, crystalseeds=crystalcourierseeds)

commutation(c6 * wanniercourier.inspace, wanniercourier' * wanniercourier) |> maximum

commutation(c6 * blockedH.outspace, blockedH) |> maximum

c6R = c6 * blockedH.outspace

c6R * blockedH * c6R' - blockedH |> maximum

c6R' * c6R - idmap(blockedH.outspace, blockedH.outspace) |> maximum

wanniercourier' * wanniercourier - idmap(wanniercourier.inspace, wanniercourier.inspace) |> maximum

(wanniercourier * wanniercourier') * wanniercourier - wanniercourier |> maximum

c6WC = c6 * wanniercourier.inspace

c6WC' * c6WC - idmap(c6WC.inspace, c6WC.inspace) |> maximum
c6WC * c6WC' - idmap(c6WC.inspace, c6WC.inspace) |> maximum

Base.:abs(fockmap::FockMap)::FockMap = FockMap(fockmap.outspace, fockmap.inspace, abs.(fockmap |> rep))

visualize(c6WC, rowrange=1:64, colrange=1:64)

c6WC * wanniercourier' - wanniercourier' |> maximum

(c6WC.outspace |> orderedmodes)[37].attrs
c6WC.outspace |> unitcellfock |> modeattrs
(c6 * (c6WC.outspace |> unitcellfock)).outspace |> modeattrs
Ft = fourier(c6WC.outspace, (c6 * (c6WC.outspace |> unitcellfock)).outspace)
visualize(FockMap(Ft.outspace, Ft.inspace, abs.(Ft |> rep)))

(wanniercourier * c6WC * wanniercourier') * wanniercourier * c6WC' - wanniercourier |> maximum

c6R' * c6R - idmap(blockedH.outspace, blockedH.outspace) |> maximum
c6R * c6R' - idmap(blockedH.outspace, blockedH.outspace) |> maximum

(c6 * wanniercourier.outspace) * wanniercourier * (c6 * wanniercourier.inspace)' - wanniercourier |> maximum

commutation(c6 * wanniercourier.inspace, wanniercourier' * wanniercourier) |> maximum

wanniercourier' * wanniercourier - idmap(wanniercourier.inspace, wanniercourier.inspace) |> maximum

courierH = wanniercourier' * blockedH * wanniercourier
courierH |> crystalspectrum |> visualize

commutation(c6 * courierH.outspace, courierH) |> maximum

fullcourierseed - (c3 * fullcourierseed.outspace) * fullcourierseed * (c3 * fullcourierseed.inspace)' |> maximum
courierseed2 - (c3 * courierseed2.outspace) * courierseed2 |> maximum
courierseedingcenter |> euclidean

visualize(courierseed |> columnspec)
visualize((c3 * courierseed.outspace) * courierseed |> columnspec)
courierseed |> makezero |> columnspec
FockMap((c3 * courierseed.outspace) * courierseed, outspace=courierseed.outspace) |> makezero |> columnspec

courierseed - (c3 * courierseed.outspace) * courierseed |> makezero |> columnspec

Ft = fourier(blockedH.outspace, blockedH.outspace |> unitcellfock)
Ft' * Ft |> visualize
relativephase(FockMap((c3 * courierseed.outspace) * courierseed, outspace=courierseed.outspace), courierseed)

modeattrs(fockspace::FockSpace)::Vector{Dict} = [m.attrs for m in fockspace |> orderedmodes]

courierseed.outspace |> unitcellfock |> modeattrs
(c3 * (courierseed.outspace |> unitcellfock)).outspace |> modeattrs
c6rep = c6 * courierseed.outspace

cbz = (courierseed.outspace |> getcrystal |> brillouinzone)
cbz == Subset(p |> basispoint for p in c3 * cbz)
visualize(cbz, Subset(p |> basispoint for p in c3 * cbz), visualspace=euclidean(RealSpace, 2))
c3c = c3 * courierseed.outspace

Base.:^(fockmap::FockMap, exponent::Integer)::FockMap = reduce(*, repeat([fockmap], exponent))
c3id = c3c^3
visualize(FockMap(c3id, outspace=c3id.inspace), rowrange=1:256, colrange=1:256)

anothercourierseed = (c6 * courierseed.outspace) * courierseed * (c6 * courierseed.inspace)'
anothercourierseedproj = anothercourierseed * anothercourierseed'
commutation(c3 * anothercourierseedproj.outspace, anothercourierseedproj) |> maximum

courierseedproj = courierseed * courierseed'
commutation(c3 * courierseedproj.outspace, courierseedproj) |> maximum

[m.attrs for m in courierseed.inspace |> orderedmodes]

courierseedproj = fullcourierseed * fullcourierseed'
commutation(c6 * courierseedproj.outspace, courierseedproj) |> maximum


function Renormalization.wannierprojection(; crystalisometries::Dict{Momentum, FockMap}, crystal::Crystal, crystalseeds::Dict{Momentum, FockMap}, svdorthothreshold::Number = 1e-1)
    wannierunitcell::Subset{Offset} = Subset(mode |> getattr(:pos) for mode in (crystalseeds |> first |> last).inspace |> orderedmodes)
    wanniercrystal::Crystal = Crystal(wannierunitcell, crystal.sizes)
    overlaps = ((k, isometry, isometry' * crystalseeds[k]) for (k, isometry) in crystalisometries)
    precarioussvdvalues::Vector = []
    function approximateisometry(k::Momentum, isometry::FockMap, overlap::FockMap)::FockMap
        U, Σ, Vt = overlap |> svd
        minsvdvalue::Number = minimum(v for (_, v) in Σ)
        if minsvdvalue < svdorthothreshold
            push!(precarioussvdvalues, minsvdvalue)
        end
        unitary::FockMap = U * Vt
        approximated::FockMap = isometry * unitary

        return FockMap(approximated, inspace=FockSpace(approximated.inspace |> orderedmodes |> setattr(:offset => k)), performpermute=false)
    end
    if (precarioussvdvalues |> length) > 0
        @warn "Precarious wannier projection with minimum svdvalue of $(precarioussvdvalues |> minimum)"
    end
    wannierisometry::FockMap = directsum(approximateisometry(k, isometry, overlap) for (k, isometry, overlap) in overlaps)
    inspace::FockSpace = FockSpace(wannierisometry.inspace, reflected=wanniercrystal)
    outspace::FockSpace = FockSpace(wannierisometry.outspace, reflected=crystal)
    return FockMap(wannierisometry, inspace=inspace, outspace=outspace, performpermute=false)
end



filledwannierisometry = wannierprojection(crystalisometries=distillresult[:filled].eigenvectors, crystal=newcrystal, seed=frozenseed)
filledwannierproj = filledwannierisometry * filledwannierisometry'
filledwannierproj |> crystalspectrum |> visualize
commutation(c6 * filledwannierproj.outspace, filledwannierproj) |> maximum
filledwanniercrystal = filledwannierisometry.inspace |> getcrystal
filledwanniercrystal.unitcell
filledbz = filledwanniercrystal |> brillouinzone
visualize(filledbz, Subset(p |> basispoint for p in c3 * filledbz), visualspace=euclidean(RealSpace, 2))
filledbz == Subset(p |> basispoint for p in c6 * filledbz)

filledH = filledwannierisometry' * blockedH * filledwannierisometry
filledH |> crystalspectrum |> visualize

c6f = (c6 * filledH.outspace)
c6filledH = c6f * filledH * c6f'
c6filledH |> crystalspectrum |> visualize

commutation(c6 * filledH.outspace, filledH) |> maximum

filledU = eigvecsh(filledH)
filledC = filledwannierisometry' * blockresult[:correlations] * filledwannierisometry
filledC |> crystalspectrum |> visualize
visualize(filledC, rowrange=1:128, colrange=1:128)

commutation(c6 * filledC.outspace, filledC) |> maximum

nakedfilledH = FockMap(filledU' * filledH * filledU, )
visualize(nakedfilledH, rowrange=1:128, colrange=1:128)
commutation(c6 * nakedfilledH.outspace, nakedfilledH) |> maximum

visualize(filledwannierisometry, rowrange=1:256, colrange=1:256)

_filledH = filledwannierproj * blockedH * filledwannierproj'
commutation(c6 * _filledH.outspace, _filledH) |> maximum

filledwannierid = filledwannierisometry' * filledwannierisometry
visualize(filledwannierid, rowrange=1:256, colrange=1:256)
commutation(c6 * filledwannierid.outspace, filledwannierid) |> maximum

c6i = c6 * filledwannierisometry.outspace
c6i' * c6i |> rep |> collect |> sum
c6f = c6 * filledwannierisometry.inspace
(c6f * filledH * c6f' - filledH) |> maximum
commutation(c6f, c6f * c6f') |> maximum
visualize(FockMap(c6f; inspace=c6f.outspace), rowrange=1:256, colrange=1:256)

commutation(c6 * blockedH.outspace, blockedH) |> maximum
commutation(c6 * filledH.outspace, filledH) |> maximum

[m.attrs for m in filledH.outspace |> unitcellfock |> orderedmodes]
visualize(fourier(filledH.inspace, filledH.outspace |> unitcellfock), rowrange=1:64)
commute(c6 * filledH.outspace, filledH; eps=1e-1)

emptywannierisometry = wannierprojection(crystalisometries=distillresult[:empty].eigenvectors, crystal=newcrystal, seed=frozenseed)
emptyH = emptywannierisometry' * blockedH * emptywannierisometry
emptyH |> crystalspectrum |> visualize

courierwannierisometry = wannierprojection(crystalisometries=distillresult[:courier].eigenvectors, crystal=newcrystal, seed=fullcourierseed)
courierH = courierwannierisometry' * blockedH * courierwannierisometry
courierH |> crystalspectrum |> visualize

commute(c6 * blockedH.outspace, blockedH; eps=5e-7)

commute(c6 * courierH.outspace, courierH; eps=1e-4)

c6H = c6 * courierH.outspace
(c6H * courierH - courierH * c6H) |> maximum

function roundingpurification(correlationspectrum::CrystalSpectrum)
    function fixeigenvalues(eigenvalue::Number)
        if isapprox(eigenvalue, 1, atol=1e-1)
            return 1
        end
        if isapprox(eigenvalue, 0, atol=1e-1)
            return 0
        end
        @warn("Eigenvalue $(eigenvalue) is not close to 0 or 1!")
    end
    eigenvalues::Dict{Mode, Number} = Dict(m => v |> fixeigenvalues for (m, v) in correlationspectrum.eigenvalues)
    return CrystalSpectrum(correlationspectrum.crystal, correlationspectrum.eigenmodes, eigenvalues, correlationspectrum.eigenvectors)
end

function crystalfockmapfromspectrum(crystalspectrum::CrystalSpectrum)::FockMap
    function momentumfockmap(k::Momentum)
        modes::Subset{Mode} = crystalspectrum.eigenmodes[k]
        eigenfock::FockSpace = modes |> FockSpace
        diagonal::FockMap = FockMap(eigenfock, eigenfock, Dict((m, m) => crystalspectrum.eigenvalues[m] |> ComplexF64 for m in modes))
        return crystalspectrum.eigenvectors[k] * diagonal * crystalspectrum.eigenvectors[k]'
    end
    fockmap::FockMap = directsum(k |> momentumfockmap for (k, _) in crystalspectrum.eigenmodes)
    crystalfock::FockSpace = FockSpace(fockmap.inspace, reflected=crystalspectrum.crystal)
    return FockMap(fockmap, inspace=crystalfock, outspace=crystalfock, performpermute=false)
end

courierC = wanniercourier' * blockresult[:correlations] * wanniercourier
purifiedcourierspectrum = courierC |> crystalspectrum |> roundingpurification
purifiedcourierspectrum |> visualize
purifiedcourierC = crystalfockmapfromspectrum(purifiedcourierspectrum)

commute(c6 * purifiedcourierC.outspace, purifiedcourierC; eps=5e-7)

block2 = blocking(:action => scale, :correlations => purifiedcourierC, :crystal => (courierC.inspace |> getcrystal))
block2H = block2[:transformer] * courierH * block2[:transformer]'
block2H |> crystalspectrum |> visualize
block2[:correlations] |> crystalspectrum |> visualize

commute(c6 * block2H.outspace, block2H; eps=5e-7)

crystalC6 = c6 * block2H.inspace

c62H = crystalC6 * block2H * crystalC6'
c62H |> crystalspectrum |> visualize

c62C = crystalC6 * block2[:correlations] * crystalC6'
c62C |> crystalspectrum |> visualize

crystal2 = block2[:crystal]

crystalpoints2::Subset{Offset} = latticepoints(crystal2)
newmodes2::Subset{Mode} = quantize(:pos, crystal2.unitcell, 1)
physicalmodes2::Subset{Mode} = spanoffset(newmodes2, crystalpoints2)
Subset(m |> pos |> euclidean for m in newmodes2) |> visualize

frozenseedingregion2::Subset{Mode} = circularregionmodes(triangular |> origin, physicalmodes2, 4.0)
Subset(m |> pos |> euclidean for m in frozenseedingregion2) |> visualize
frozenseedingfock2::FockSpace = FockSpace(frozenseedingregion2)

regionC6 = c6 * frozenseedingregion2

localC = Renormalization.regioncorrelations(block2[:correlations], frozenseedingfock2)

(regionC6 * localC * regionC6' - localC)

visualize(c6 * frozenseedingregion2)

isometries = localfrozenisometries(block2[:correlations], frozenseedingfock2, selectionstrategy=frozenselectionbycount(3))
localfilledproj = isometries[:filled] * isometries[:filled]'
C6localfilledproj = regionC6 * localfilledproj * regionC6'
FockMap(C6localfilledproj, inspace=localfilledproj.inspace, outspace=localfilledproj.outspace) |> visualize
visualize(localfilledproj)
visualize(regionC6 * localfilledproj * regionC6')

region2 = Subset(m |> pos for m in frozenseedingregion2)
c6 * region2 == region2

visualize(block2[:correlations], colrange=1:128, rowrange=1:128)

function _globaldistillerhamiltonian(;
    correlations::FockMap, restrictspace::FockSpace, localisometryselectionstrategy, manualeigenenergies::Dict{Symbol, <:Number} = Dict(:filled => -1, :empty => 1),
    symmetry::AffineTransform)

    localisometries::Dict{Symbol} = localfrozenisometries(correlations, restrictspace, selectionstrategy=localisometryselectionstrategy)
    crystalprojectors::Dict{Symbol, FockMap} = Dict(
        name => crystalprojector(localisometry=localisometries[name], crystalfock=correlations.inspace)
        for (name, isometry) in localisometries)
    globaldistillhamiltonian::FockMap = reduce(+, manualeigenenergies[name] * crystalprojector for (name, crystalprojector) in crystalprojectors)

    transformers::Vector{FockMap} = [transformation * globaldistillhamiltonian.outspace for transformation in symmetry |> pointgroupelements]

    return reduce(+, (transformer * globaldistillhamiltonian * transformer' for transformer in transformers))
end

localFt = fourier(block2[:correlations].inspace, frozenseedingfock2) / sqrt(256)
localC = localFt' * block2[:correlations] * localFt
c6localt = c6 * localC.outspace

visualize(localC)
visualize(c6localt * localC * c6localt')

localcvals = eigvalsh(localC)
plot(scatter(y=map(p -> p.second, localcvals), mode="markers"))
localiso2 = localfrozenisometries(block2[:correlations], frozenseedingfock2, selectionstrategy=frozenselectionbycount(5))
globalfilled2 = crystalprojector(localisometry=localiso2[:filled], crystalfock=block2[:correlations].inspace)
filledproj2 = globalfilled2 * globalfilled2'
filledproj2 |> crystalspectrum |> visualize

c6t = c6 * filledproj2.outspace
(c6t * filledproj2 * c6t') |> crystalspectrum |> visualize

dH2 = globaldistillerhamiltonian(;
    correlations=block2[:correlations],
    restrictspace=frozenseedingfock2,
    localisometryselectionstrategy=frozenselectionbycount(3),
    symmetry=c6)

dH2 |> crystalspectrum |> visualize
plot([v for (_, v) in  ((c6t * dH2 * c6t') |> crystalspectrum).eigenvalues] |> sort)
(c6 * (dH2.outspace |> unitcellfock)) |> visualize

Ft = fourier(dH2.outspace, frozenseedingfock2) / sqrt(256)
visualize(Ft, rowrange=1:256)

visualize(Ft' * block2[:correlations] * Ft)
corrvals2, _ = eigh(Ft' * block2[:correlations] * Ft)
plot([v for (_, v) in corrvals2])

c6t = c6 * dH2.outspace
visualize(c6t * c6t, rowrange=1:256, colrange=1:256)
c66 = c6t * c6t * c6t * c6t * c6t * c6t
visualize(FockMap(c66, inspace=c66.outspace), rowrange=1:256, colrange=1:256)
dH2 |> crystalspectrum |> visualize
(c6t * dH2 * c6t') |> crystalspectrum |> visualize

[m.attrs for m in c6t.outspace |> unitcellfock |> orderedmodes]

distillresult2 = dH2 |> distillation

bigregion::Subset{Mode} = filter(circularfilter(origin(euclidean(RealSpace, 2)), 2.0), physicalmodes)
bigregionfock::FockSpace = FockSpace(bigregion)
Ft = fourier(blockedH.outspace, bigregionfock)

filledisometries = Ft' * filledwannierisometry
columns(filledisometries, (filledisometries.inspace |> orderedmodes)[1] |> Subset |> FockSpace) |> columnspec |> visualize

overlaps = [unitary' * rows(seed, unitary.outspace) for (k, unitary) in distillresult[:filled].eigenvectors]
u, s, vt = overlaps[1] |> svd
[s...]

visualize(u * vt)

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

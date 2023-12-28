using Plotly, SmithNormalForm, LinearAlgebra, OrderedCollections, SparseArrays, Combinatorics
using Revise
using Zipper
# include("../src/Zipper.jl")
# using .Zipper


# function findeigenfunction(transformation::AffineTransform; eigenvalue::Number = 1)::BasisFunction
#     eigenvaluekey::Tuple = hashablecomplex(eigenvalue |> Complex, transformation.eigenvaluehashdenominator)
#     if !haskey(transformation.eigenfunctions, eigenvaluekey)
#         error("No eigenfunction is found for eigenvalue $(eigenvalue)!")
#     end
#     return transformation.eigenfunctions[eigenvaluekey]
# end

triangular = RealSpace([sqrt(3)/2 -1/2; 0. 1.]')
kspace = convert(MomentumSpace, triangular)

pa = [1/3, 2/3] ∈ triangular
pb = [2/3, 1/3] ∈ triangular
pc = (pa + pb) / 2
spatialsnappingcalibration((pa, pb, pc))

c6 = pointgrouptransform([cos(π/3) -sin(π/3); sin(π/3) cos(π/3)])

unitcell = Subset(pa, pb)
crystal = Crystal(unitcell, [96, 96])
reciprocalhashcalibration(crystal.sizes)

# test for small lattice##########################################################################################
L=6
crystaltest = Crystal(unitcell, [L, L])
energyspectrumtest = computeenergyspectrum(bonds, crystal=crystaltest)
# energyspectrumtest |> visualize
groundstatestest = groundstatespectrum(energyspectrumtest, perunitcellfillings=1)
# groundstatestest |> visualize
Chern_number_multiband(groundstatestest)


groundstateprojectortest = groundstatestest |> crystalprojector
Ctest = idmap(groundstateprojectortest.outspace) - groundstateprojectortest
Ctest |> visualize

new_C_test = blocked_correlation(Ctest,1) # 2^n
new_C_test |>visualize

new_gs_test = groundstatespectrum(new_C_test |> crystalspectrum, perunitcellfillings = 4)
# new_gs_test|> geteigenvectors
# new_gs_test |> visualize
Chern_number_multiband(new_gs_test)


# understand the current blocking function
function Base.:*(scale::Scale, crystalfock::FockSpace{Crystal})::FockMap
    crystal::Crystal = crystalfock |> getcrystal
    scaledcrystal::Crystal = scale * crystal
    unscaledblockedregion::Subset{Offset} = (scale |> inv) * scaledcrystal.unitcell
    bz::Subset{Momentum} = crystal |> brillouinzone
    basismodes::Subset{Mode} = crystalfock |> rep |> first
    scaledbz::Subset{Momentum} = scaledcrystal |> brillouinzone

    momentummappings::Base.Generator = (basispoint(scale * p) => p for p in bz)     #?
    mappingpartitions::Dict{Point, Vector{Point}} = foldl(momentummappings; init=Dict{Point, Vector{Point}}()) do d, (k, v)
        mergewith!(append!, d, LittleDict(k => [v]))
    end         #?

    ksubsets::Dict{Momentum, Subset{Mode}} = crystalfock |> crystalsubsets      # c(k) dictionary
    scaledfock::FockSpace = ((ksubsets[k] for k in mappingpartitions[kscaled]) |> subsetunion |> FockSpace
        for kscaled in scaledbz) |> fockspaceunion          # scaled Fock space

    # TODO: Use of latticeoff requires revisit, since latticeoff only truncates the decimals.
    unscaledblockedlatticeoffsets::Subset{Offset} = Subset(pbc(crystal, p) |> latticeoff for p in unscaledblockedregion)        # 
    unscaledblockedunitcellfock::FockSpace = spanoffset(basismodes, unscaledblockedlatticeoffsets)

    restrictedfourier::FockMap = fourier(bz, unscaledblockedunitcellfock)     # fourier(momentum points, FockSpace) gives a Fock Map representing a Fourier transform from the given FockSpace
    volumeratio::Number = (crystal |> vol) / (scaledcrystal |> vol)
    permutedfourier::FockMap = Zipper.permute(restrictedfourier, outspace=scaledfock) / sqrt(volumeratio)  # Rescale the FT due to space rescaling

    function repackfourierblocks(source::FockMap, kscaled::Momentum, partition::Subset{Mode})::FockMap
        partitionrows::FockMap = rows(source, partition |> FockSpace)
        inspace::FockSpace = (Subset(setattr(mode, :offset => kscaled, :pos => scale * convert(Point, mode))
            for mode in partitionrows.inspace |> orderedmodes)
            |> FockSpace)
        return FockMap(partitionrows.outspace, inspace, partitionrows |> rep)
    end

    repackedblocks::Base.Generator = (
        repackfourierblocks(permutedfourier, kscaled, partition)
        for (kscaled, partition) in Iterators.zip(scaledbz, scaledfock |> rep))
    blocking::FockMap = directsum(repackedblocks)
    return FockMap(blocking, inspace=FockSpace(blocking |> getinspace, reflected=scaledcrystal), outspace=crystalfock)'
end

function correlation_blocking(correlations::FockMap, scale)
    initial_crystal = correlations |> getoutspace |> getcrystal
    scaled_crystal::Crystal = scale * initial_crystal
    unscaledblockedregion::Subset{Offset} = (scale |> inv) * scaledcrystal.unitcell         
    bz::Subset{Momentum} = crystal |> brillouinzone
    basismodes::Subset{Mode} = crystalfock |> rep |> first
    scaledbz::Subset{Momentum} = scaledcrystal |> brillouinzone

    momentummappings::Base.Generator = (basispoint(scale * p) => p for p in bz)     #?
    mappingpartitions::Dict{Point, Vector{Point}} = foldl(momentummappings; init=Dict{Point, Vector{Point}}()) do d, (k, v)
        mergewith!(append!, d, LittleDict(k => [v]))
    end         #?
end
##############################################################################################################################

modes::Subset{Mode} = quantize(:pos, unitcell, 1)
m0, m1 = members(modes)
# sparsefock(modes, crystal|> latticepoints)
tₙ = -1 + 0im
tₕ = 0.1im

nearestneighbor = [
    (m0, m1) => tₙ,
    (m0, m1 |> setattr(:offset => [-1, 0] ∈ triangular)) => tₙ,
    (m0, m1 |> setattr(:offset => [0, 1] ∈ triangular)) => tₙ]

haldane = [
    (m0, m0 |> setattr(:offset => [1, 1] ∈ triangular)) => tₕ,
    (m0, m0 |> setattr(:offset => [-1, 0] ∈ triangular)) => tₕ,
    (m0, m0 |> setattr(:offset => [0, -1] ∈ triangular)) => tₕ,
    (m1, m1 |> setattr(:offset => [1, 1] ∈ triangular)) => -tₕ,
    (m1, m1 |> setattr(:offset => [-1, 0] ∈ triangular)) => -tₕ,
    (m1, m1 |> setattr(:offset => [0, -1] ∈ triangular)) => -tₕ]

bonds::FockMap = bondmap([nearestneighbor..., haldane...])

energyspectrum = computeenergyspectrum(bonds, crystal=crystal)
energyspectrum |> visualize

groundstates::CrystalSpectrum = groundstatespectrum(energyspectrum, perunitcellfillings=1)
groundstates |> visualize

groundstateprojector = groundstates |> crystalprojector
# correlation = 1 - projector
C = idmap(groundstateprojector.outspace) - groundstateprojector

C |> crystalspectrum |> visualize

inplaceadd(a::FockMap, b::FockMap)::FockMap = a + FockMap(b, outspace=a|>getoutspace, inspace=a|>getinspace, performpermute=false)

entanglemententropy(C |> crystalspectrum) 

# groundstatenew = groundstatespectrum(C|> crystalspectrum, perunitcellfillings=1.0)

statevectors::Dict{Momentum, FockMap} = groundstates|>geteigenvectors

function computequantumgeometrictensor2d(state::CrystalSpectrum)
    statevectors::Dict{Momentum, FockMap} = state|>geteigenvectors
    statecrystal = state |> getcrystal
    Δkx = (1 / size(statecrystal)[1], 0) ∈ kspace
    Δky = (0, 1 / size(statecrystal)[2]) ∈ kspace

    kprojectors::Dict{Momentum, FockMap} = Dict(k => kstate * kstate' for (k, kstate) in statevectors)
    kcorrelations::Dict{Momentum, FockMap} = Dict(k => idmap(projector|>getoutspace) - projector for (k, projector) in kprojectors)  # correlations |> crystalsubmaps |> Dict
    ΔCkxs::Dict{Momentum, FockMap} = Dict(k => inplaceadd(-kcorr, kcorrelations[k + Δkx |> basispoint]) for (k, kcorr) in kcorrelations)
    ΔCkys::Dict{Momentum, FockMap} = Dict(k => inplaceadd(-kcorr, kcorrelations[k + Δky |> basispoint]) for (k, kcorr) in kcorrelations)
    
    function kqgtensor(k::Momentum)::FockMap
        gxx = statevectors[k]' * (ΔCkxs[k]' * kcorrelations[k] * ΔCkxs[k]) * statevectors[k]
        gxy = statevectors[k]' * (ΔCkxs[k]' * kcorrelations[k] * ΔCkys[k]) * statevectors[k]
        gyx = statevectors[k]' * (ΔCkys[k]' * kcorrelations[k] * ΔCkxs[k]) * statevectors[k]
        gyy = statevectors[k]' * (ΔCkys[k]' * kcorrelations[k] * ΔCkys[k]) * statevectors[k]

        xmode::Mode = Mode(:offset=>k, :axis=>(:x))
        ymode::Mode = Mode(:offset=>k, :axis=>(:y))

        fockspace::FockSpace = FockSpace([xmode, ymode])
        return FockMap(fockspace, fockspace, [rep(gxx)[1,1] rep(gxy)[1,1]; rep(gyx)[1,1] rep(gyy)[1,1]])
    end

    return Dict(k=>k|>kqgtensor for (k, _) in statevectors)
end

Base.:real(fockmap::FockMap)::FockMap = FockMap(fockmap|>getoutspace, fockmap|>getinspace, fockmap|>rep|>real)
Base.:imag(fockmap::FockMap)::FockMap = FockMap(fockmap|>getoutspace, fockmap|>getinspace, fockmap|>rep|>imag)

# obsolete
function get_Chern_number(state::CrystalSpectrum)
    kqgtensors = computequantumgeometrictensor2d(state)
    gxx = directsum(tensor[1,2] for (k, tensor) in kqgtensors)
    barecrystal = Crystal(triangular|>getorigin, crystal|>size)
    barecrystalfock = FockSpace(gxx|>getoutspace, reflected=barecrystal)
    bc = directsum(tensor[1,2]|>imag for (k, tensor) in kqgtensors)
    barecrystalfock = FockSpace(bc|>getoutspace, reflected=barecrystal)
    berrycurvature = 2 * FockMap(bc, inspace=barecrystalfock, outspace=barecrystalfock, performpermute=false)

    return (berrycurvature|>rep|>sum) / (2π)
end

get_Chern_number(groundstates)

# https://arxiv.org/pdf/cond-mat/0503172.pdf
function Chern_number(state::CrystalSpectrum)
    statevectors::Dict{Momentum, FockMap} = state|>geteigenvectors
    statecrystal = state |> getcrystal
    Δkx = (1 / size(statecrystal)[1], 0) ∈ kspace
    Δky = (0, 1 / size(statecrystal)[2]) ∈ kspace
    Δk = [Δkx, Δky]
    
    function LinkVariable(k::Momentum, direction)
        delta_k = Δk[direction]
        U_mu = (statevectors[k] |> rep)' * (statevectors[k+delta_k|> basispoint]|> rep)
        U_mu = U_mu[1]
        U_mu /= abs(U_mu)
        return U_mu
    end
    
    function Berry_phase(k::Momentum)
        return LinkVariable(k, 1) * LinkVariable(k + Δk[1]|> basispoint, 2) / (LinkVariable(k + Δk[2]|> basispoint, 1) * LinkVariable(k, 2))
    end

    function Berry_curvature(k::Momentum)
        return -1im*log(Berry_phase(k))
    end

    C = 0 
    for kpoints in statevectors
        k = kpoints[1]
        C += Berry_curvature(k)
    end
    return C/(2*pi)
end

function Chern_number_multiband(state::CrystalSpectrum)
    statevectors::Dict{Momentum, FockMap} = state |> geteigenvectors
    statecrystal = state |> getcrystal
    Δkx = (1 / size(statecrystal)[1], 0) ∈ kspace
    Δky = (0, 1 / size(statecrystal)[2]) ∈ kspace
    Δk = [Δkx, Δky]
    
    function LinkVariable(k::Momentum, direction)
        delta_k = Δk[direction]
        println("statevector value k: ", statevectors[k] |> rep)
        println("statevector value k': ",statevectors[k+delta_k|> basispoint] |> rep)
        overlap = (statevectors[k] |> rep)' * (statevectors[k+delta_k|> basispoint]|> rep)
        U, Σ, Vt = Matrix(overlap) |> svd
        println("sigma mat: ",Σ)
        U_mu = det(U*Vt)
        U_mu /= abs(U_mu)
        println("U val: ",U_mu)
        return U_mu
    end
    
    function Berry_phase(k::Momentum)
        return LinkVariable(k, 1) * LinkVariable(k + Δk[1]|> basispoint, 2) / (LinkVariable(k + Δk[2]|> basispoint, 1) * LinkVariable(k, 2))
    end

    function Berry_curvature(k::Momentum)
        return -1im*log(Berry_phase(k))
    end

    C = 0 
    for kpoints in statevectors
        k = kpoints[1]
        C += Berry_curvature(k)
    end
    return C/(2*pi)
end
Chern_number_multiband(new_gs_test)
Chern_number_multiband(groundstatestest)


# function Berry_curvature(state::CrystalSpectrum)
#     statevectors::Dict{Momentum, FockMap} = state|>geteigenvectors
#     statecrystal = state |> getcrystal
#     Δkx = (1 / size(statecrystal)[1], 0) ∈ kspace
#     Δky = (0, 1 / size(statecrystal)[2]) ∈ kspace
#     Δk = [Δkx, Δky]
    
#     function LinkVariable(k::Momentum, direction)
#         delta_k = Δk[direction]
#         println("kin", statevectors[k]' |> getoutspace)
#         # println("kout", statevectors[k] |> getoutspace)
#         println("kpin",statevectors[k+delta_k] |> getinspace)
#         # println("kpout",statevectors[k+delta_k] |> getoutspace)
#         Ikk = idmap(statevectors[k] |> getinspace, statevectors[k+delta_k] |> getinspace)
#         # U_mu = (statevectors[k] |> rep)' * (statevectors[k+delta_k|> basispoint]|> rep)
#         U_mu = statevectors[k]' * statevectors[k+delta_k|> basispoint] 
#         # U_mu = U_mu[1]
#         # print(U_mu)
#         norm = abs(rep(U_mu)[1])
#         U_mu /= norm
#         return U_mu
#     end

#     # function LinkVariable(k::Momentum, direction)
#     #     delta_k = Δk[direction]
#     #     U_mu = (statevectors[k] |> rep)' * (statevectors[k+delta_k|> basispoint]|> rep)
#     #     U_mu = U_mu[1]
#     #     U_mu /= abs(U_mu)
#     #     return U_mu
#     # end
    
#     function Berry_phase(k::Momentum)
#         return LinkVariable(k, 1) * LinkVariable(k + Δk[1]|> basispoint, 2) / (LinkVariable(k + Δk[2]|> basispoint, 1) * LinkVariable(k, 2))
#     end
#     return Dict(k |> vec => real(-1im*log(Berry_phase(k))) for (k,_) in statevectors)
#     # bc = directsum(-1im*log(Berry_phase(k)) for (k, _) in statevectors)
#     # barecrystal = Crystal(triangular|>getorigin, statecrystal|>size)
#     # barecrystalfock = FockSpace(bc|>getoutspace, reflected=barecrystal)
#     # return FockMap(bc, inspace=barecrystalfock, outspace=barecrystalfock, performpermute=false)
# end

# function Chern_number(state::CrystalSpectrum)
#     statevectors::Dict{Momentum, FockMap} = state|>geteigenvectors
#     C = 0 
#     for kpoints in statevectors
#         k = kpoints[1]
#         print(k)
#         C += BerryCurvature(state)[:k]
#     end
#     return C 
# end
# Berry_cur = Berry_curvature(groundstates)

# x_points = [x for (x,y) in keys(Berry_cur)]
# y_points = [y for (x,y) in keys(Berry_cur)]
# z_points = [bc for (k, bc) in Berry_cur]

# heatmap(x_points, y_points, z_points)

# Berry_curvature(groundstates)
Chern_number(groundstates)

statevectors::Dict{Momentum, FockMap} = groundstates|>geteigenvectors
statevectors

function zipper0(correlations::FockMap)
    crystalfock = correlations.outspace

    scale = Scale([2 0; 0 2], crystalfock |> getcrystal |> getspace)
    blockresult = blocking(:action => scale, :correlations => correlations, :crystal => crystalfock |> getcrystal)

    blockedcrystal::Crystal = blockresult[:crystal]
    blockedcorrelations::FockMap = blockresult[:correlations]
    blocker = blockresult[:transformer]

    function circularregionmodes(origin::Offset, physicalmodes::Subset{Mode}, radius::Number)::Subset{Mode}
        currentspace::RealSpace = correlations.outspace |> getcrystal |> getspace |> orthogonalspace
        physicalnorm = m -> lineartransform(currentspace, m |> getpos) |> norm
        return filter(p -> physicalnorm(p - origin) < radius, physicalmodes)
    end

    crystalpoints::Subset{Offset} = latticepoints(blockedcrystal)
    samplepoints::Subset{Offset} = crystalpoints + c6^2 * crystalpoints + c6^4 * crystalpoints
    blockedmodes::Subset{Mode} = quantize(:pos, blockedcrystal.unitcell, 1)
    physicalmodes::Subset{Mode} = spanoffset(blockedmodes, samplepoints)


    frozenseedingmodes::Subset{Mode} = circularregionmodes(triangular |> getorigin, physicalmodes, 2)
    # visualize(frozenseedingregion, title="Frozen Seeding Region", visualspace=euclidean(RealSpace, 2))
    frozenseedingfock::FockSpace = FockSpace{Region}(frozenseedingmodes)

    globaldistiller = globaldistillerhamiltonian(
        correlations=blockresult[:correlations],
        restrictspace=frozenseedingfock,
        localisometryselectionstrategy=frozenselectionbycount(3))

    globaldistillerspectrum = globaldistiller |> crystalspectrum
    # visualize(globaldistillerspectrum, title="Global Distiller")

    distillresult = distillation(globaldistillerspectrum, :courier => v -> abs(v) < 1e-5, :empty => v -> v > 1e-5, :filled => v -> v < -1e-5)

    courierseedingcenter::Offset = [2/3, 1/3] ∈ (blockedmodes |> getspace)
    courierseedingmodes::Subset{Mode} = circularregionmodes(courierseedingcenter, physicalmodes, 1.8)
    courierseedingregion::Subset{Offset} = Subset(m |> getpos for m in courierseedingmodes)
    # visualize(courierseedingregion, courierseedingcenter |> Subset, title="Courier Seeding Region", visualspace=euclidean(RealSpace, 2))
    courierseedingfock::FockSpace{Region} = FockSpace{Region}(courierseedingmodes)

    c3 = c6^2 |> recenter(courierseedingcenter)

    blockedcourierprojector = distillresult[:courier] |> crystalprojector
    blockedcouriercorrelation = idmap(blockedcourierprojector.outspace, blockedcourierprojector.outspace) - blockedcourierprojector

    localcourierseed = findlocalspstates(statecorrelations=blockedcouriercorrelation, regionfock=courierseedingfock, symmetry=c3, spectrumextractpredicate=v -> v < 5e-2)[1]
    fullcourierseed = localcourierseed + (c6 * localcourierseed.outspace) * localcourierseed * (c6 * localcourierseed.inspace)'

    crystalcourierseeds = crystalisometries(localisometry=fullcourierseed, crystalfock=blockedcorrelations.outspace, addinspacemomentuminfo=true)

    wanniercourierisometry = wannierprojection(
        crystalisometries=distillresult[:courier].eigenvectors, crystal=blockedcrystal, crystalseeds=crystalcourierseeds)

    # regionalrestriction(wanniercourierisometry, courierseedingfock) |> visualize

    couriercorrelations = wanniercourierisometry' * blockedcorrelations * wanniercourierisometry
    # couriercorrelations |> crystalspectrum |> visualize
    couriercorrelationspectrum = couriercorrelations |> crystalspectrum
    purifiedcorrelationspectrum = couriercorrelationspectrum |> roundingpurification
    # purifiedcorrelationspectrum |> visualize
    unpurecorrelations = couriercorrelations
    couriercorrelations = purifiedcorrelationspectrum |> FockMap

    blockedfilledprojector = distillresult[:filled] |> crystalprojector
    blockedfilledcorrelation = idmap(blockedfilledprojector.outspace, blockedfilledprojector.outspace) - blockedfilledprojector
    filledseed = findlocalspstates(statecorrelations=blockedfilledcorrelation, regionfock=frozenseedingfock, symmetry=c6, spectrumextractpredicate=v -> v < 1e-2, degeneracythreshold=1e-3)[3]

    crystalfilledseeds = crystalisometries(localisometry=filledseed, crystalfock=blockedcorrelations.outspace, addinspacemomentuminfo=true)
    wannierfilledisometry = wannierprojection(
        crystalisometries=distillresult[:filled].eigenvectors, crystal=blockedcrystal, crystalseeds=crystalfilledseeds)

    blockedemptyprojector = distillresult[:empty] |> crystalprojector
    blockedemptycorrelation = idmap(blockedemptyprojector.outspace, blockedemptyprojector.outspace) - blockedemptyprojector
    emptyseed = findlocalspstates(statecorrelations=blockedemptycorrelation, regionfock=frozenseedingfock, symmetry=c6, spectrumextractpredicate=v -> v < 1e-2, degeneracythreshold=1e-3)[3]

    crystalemptyseeds = crystalisometries(localisometry=emptyseed, crystalfock=blockedcorrelations.outspace, addinspacemomentuminfo=true)
    wannieremptyisometry = wannierprojection(
        crystalisometries=distillresult[:empty].eigenvectors, crystal=blockedcrystal, crystalseeds=crystalemptyseeds)

    filledcorrelations = wannierfilledisometry' * blockedcorrelations * wannierfilledisometry
    emptycorrelations = wannieremptyisometry' * blockedcorrelations * wannieremptyisometry

    isometry = wannierfilledisometry + wanniercourierisometry + wannieremptyisometry    # the total isometry (unitary)
    u_zipper = blocker' * isometry

    C_total = u_zipper' * correlations * u_zipper # =filledcorrelations + couriercorrelations + emptycorrelations      # total correlation matrix of the resulting lattice

    return Dict(
        :blocker => blocker,
        :blockedcorrelations => blockedcorrelations,
        :correlations => couriercorrelations, 
        :courierzipper => wanniercourierisometry' * blocker, 
        :filledzipper => wannierfilledisometry' * blocker,
        :emptyzipper => wannieremptyisometry' * blocker,
        :globaldistiller => globaldistiller, 
        :filledcorrelations => filledcorrelations, 
        :emptycorrelations => emptycorrelations,
        :entanglemententropy => entanglemententropy(couriercorrelationspectrum),
        :wanniercourierisometry => wanniercourierisometry,
        :wannierfilledisometry => wannierfilledisometry,
        :wannieremptyisometry => wannieremptyisometry, 
        :C_total => C_total,
        :isometry => isometry,
        :u_zipper => u_zipper,
        :frozenseedingfock => frozenseedingfock,
        :courierseedingfock => courierseedingfock,
        :distillresult => distillresult)
end

function zipper(correlations::FockMap)
    crystalfock = correlations.outspace

    scale = Scale([2 0; 0 2], crystalfock |> getcrystal |> getspace)
    blockresult = blocking(:action => scale, :correlations => correlations, :crystal => crystalfock |> getcrystal)

    blockedcrystal::Crystal = blockresult[:crystal]
    blockedcorrelations::FockMap = blockresult[:correlations]
    blocker = blockresult[:transformer]

    function circularregionmodes(origin::Offset, physicalmodes::Subset{Mode}, radius::Number)::Subset{Mode}
        currentspace::RealSpace = correlations.outspace |> getcrystal |> getspace |> orthogonalspace
        physicalnorm = m -> lineartransform(currentspace, m |> getpos) |> norm
        return filter(p -> physicalnorm(p - origin) < radius, physicalmodes)
    end

    crystalpoints::Subset{Offset} = latticepoints(blockedcrystal)
    samplepoints::Subset{Offset} = crystalpoints + c6^2 * crystalpoints + c6^4 * crystalpoints
    blockedmodes::Subset{Mode} = quantize(:pos, blockedcrystal.unitcell, 1)
    physicalmodes::Subset{Mode} = spanoffset(blockedmodes, samplepoints)


    frozenseedingmodes::Subset{Mode} = circularregionmodes(triangular |> getorigin, physicalmodes, 2)
    # visualize(frozenseedingregion, title="Frozen Seeding Region", visualspace=euclidean(RealSpace, 2))
    frozenseedingfock::FockSpace = FockSpace{Region}(frozenseedingmodes)

    globaldistiller = globaldistillerhamiltonian(
        correlations=blockresult[:correlations],
        restrictspace=frozenseedingfock,
        localisometryselectionstrategy=frozenselectionbycount(3))

    globaldistillerspectrum = globaldistiller |> crystalspectrum
    # visualize(globaldistillerspectrum, title="Global Distiller")

    distillresult = distillation(globaldistillerspectrum, :courier => v -> abs(v) < 1e-5, :empty => v -> v > 1e-5, :filled => v -> v < -1e-5)

    courierseedingcenter::Offset = [2/3, 1/3] ∈ (blockedmodes |> getspace)
    courierseedingmodes::Subset{Mode} = circularregionmodes(courierseedingcenter, physicalmodes, 1.8)
    courierseedingregion::Subset{Offset} = Subset(m |> getpos for m in courierseedingmodes)
    # visualize(courierseedingregion, courierseedingcenter |> Subset, title="Courier Seeding Region", visualspace=euclidean(RealSpace, 2))
    courierseedingfock::FockSpace{Region} = FockSpace{Region}(courierseedingmodes)

    c3 = c6^2 |> recenter(courierseedingcenter)

    blockedcourierprojector = distillresult[:courier] |> crystalprojector
    blockedcouriercorrelation = idmap(blockedcourierprojector.outspace, blockedcourierprojector.outspace) - blockedcourierprojector

    localcourierseed = findlocalspstates(statecorrelations=blockedcouriercorrelation, regionfock=courierseedingfock, symmetry=c3, spectrumextractpredicate=v -> v < 5e-2)[1]
    fullcourierseed = localcourierseed + (c6 * localcourierseed.outspace) * localcourierseed * (c6 * localcourierseed.inspace)'

    crystalcourierseeds = crystalisometries(localisometry=fullcourierseed, crystalfock=blockedcorrelations.outspace, addinspacemomentuminfo=true)

    wanniercourierisometry = wannierprojection(
        crystalisometries=distillresult[:courier].eigenvectors, crystal=blockedcrystal, crystalseeds=crystalcourierseeds)

    # regionalrestriction(wanniercourierisometry, courierseedingfock) |> visualize

    couriercorrelations = wanniercourierisometry' * blockedcorrelations * wanniercourierisometry
    # couriercorrelations |> crystalspectrum |> visualize
    couriercorrelationspectrum = couriercorrelations |> crystalspectrum
    purifiedcorrelationspectrum = couriercorrelationspectrum |> roundingpurification
    # purifiedcorrelationspectrum |> visualize
    unpurecorrelations = couriercorrelations
    couriercorrelations = purifiedcorrelationspectrum |> FockMap

    blockedfilledprojector = distillresult[:filled] |> crystalprojector
    blockedfilledcorrelation = idmap(blockedfilledprojector.outspace, blockedfilledprojector.outspace) - blockedfilledprojector
    filledseed = findlocalspstates(statecorrelations=blockedfilledcorrelation, regionfock=frozenseedingfock, symmetry=c6, spectrumextractpredicate=v -> v < 1e-2, degeneracythreshold=1e-3)[3]

    crystalfilledseeds = crystalisometries(localisometry=filledseed, crystalfock=blockedcorrelations.outspace, addinspacemomentuminfo=true)
    wannierfilledisometry = wannierprojection(
        crystalisometries=distillresult[:filled].eigenvectors, crystal=blockedcrystal, crystalseeds=crystalfilledseeds)

    blockedemptyprojector = distillresult[:empty] |> crystalprojector
    blockedemptycorrelation = idmap(blockedemptyprojector.outspace, blockedemptyprojector.outspace) - blockedemptyprojector
    emptyseed = findlocalspstates(statecorrelations=blockedemptycorrelation, regionfock=frozenseedingfock, symmetry=c6, spectrumextractpredicate=v -> v < 1e-2, degeneracythreshold=1e-3)[3]

    crystalemptyseeds = crystalisometries(localisometry=emptyseed, crystalfock=blockedcorrelations.outspace, addinspacemomentuminfo=true)
    wannieremptyisometry = wannierprojection(
        crystalisometries=distillresult[:empty].eigenvectors, crystal=blockedcrystal, crystalseeds=crystalemptyseeds)

    filledcorrelations = wannierfilledisometry' * blockedcorrelations * wannierfilledisometry
    emptycorrelations = wannieremptyisometry' * blockedcorrelations * wannieremptyisometry

    isometry = wannierfilledisometry + wanniercourierisometry + wannieremptyisometry    # the total isometry (unitary)
    u_zipper = blocker' * isometry

    C_total = u_zipper' * correlations * u_zipper

    #modifiedTotalCorrelationMat = modified_filledcorrelations + correlations + modified_emptycorrelations
    
    swavefockspace = FockSpace(mode for mode in filledcorrelations |> getinspace if swave == (mode |> getorbital(swave)))
    swavefockspace = FockSpace(swavefockspace, reflected=wannierfilledisometry |> getinspace |> getcrystal)
    swavefilledisometry = wannierfilledisometry[:, swavefockspace]
    
    swavefilledprojector = swavefilledisometry * swavefilledisometry'
    swavefilledcorrelations = idmap(swavefilledprojector |> getinspace) - swavefilledprojector
    filledcorrelations_tilde = wannierfilledisometry' * swavefilledcorrelations * wannierfilledisometry
    
    C_tilde = filledcorrelations_tilde + couriercorrelations + emptycorrelations 

    C_previous_tilde = u_zipper * C_tilde * u_zipper'
    C_previous= u_zipper * C_total * u_zipper'

    return Dict(
        :blocker => blocker,
        :blockedcorrelations => blockedcorrelations,
        :correlations => couriercorrelations, 
        :courierzipper => wanniercourierisometry' * blocker, 
        :filledzipper => wannierfilledisometry' * blocker,
        :emptyzipper => wannieremptyisometry' * blocker,
        :globaldistiller => globaldistiller, 
        :filledcorrelations => filledcorrelations, 
        :emptycorrelations => emptycorrelations,
        :entanglemententropy => entanglemententropy(couriercorrelationspectrum),
        :wanniercourierisometry => wanniercourierisometry,
        :wannierfilledisometry => wannierfilledisometry,
        :wannieremptyisometry => wannieremptyisometry, 
        :C_total => C_total,
        :filledcorrelations_tilde => filledcorrelations_tilde,
        :C_tilde => C_tilde,
        :isometry => isometry,
        :u_zipper => u_zipper,
        :C_previous_tilde => C_previous_tilde,
        :C_previous => C_previous,
        :frozenseedingfock => frozenseedingfock,
        :courierseedingfock => courierseedingfock,
        :distillresult => distillresult)
end 

function zipper_even(correlations::FockMap)
    crystalfock = correlations.outspace

    scale = Scale([2 0; 0 2], crystalfock |> getcrystal |> getspace)
    blockresult = blocking(:action => scale, :correlations => correlations, :crystal => crystalfock |> getcrystal)

    blockedcrystal::Crystal = blockresult[:crystal]
    blockedcorrelations::FockMap = blockresult[:correlations]
    blocker = blockresult[:transformer]

    function circularregionmodes(origin::Offset, physicalmodes::Subset{Mode}, radius::Number)::Subset{Mode}
        currentspace::RealSpace = correlations.outspace |> getcrystal |> getspace |> orthogonalspace
        physicalnorm = m -> lineartransform(currentspace, m |> getpos) |> norm
        return filter(p -> physicalnorm(p - origin) < radius, physicalmodes)
    end

    crystalpoints::Subset{Offset} = latticepoints(blockedcrystal)
    samplepoints::Subset{Offset} = crystalpoints + c6^2 * crystalpoints + c6^4 * crystalpoints
    blockedmodes::Subset{Mode} = quantize(:pos, blockedcrystal.unitcell, 1)
    physicalmodes::Subset{Mode} = spanoffset(blockedmodes, samplepoints)


    frozenseedingmodes::Subset{Mode} = circularregionmodes(triangular |> getorigin, physicalmodes, 2)
    # visualize(frozenseedingregion, title="Frozen Seeding Region", visualspace=euclidean(RealSpace, 2))
    frozenseedingfock::FockSpace = FockSpace{Region}(frozenseedingmodes)

    globaldistiller = globaldistillerhamiltonian(
        correlations=blockresult[:correlations],
        restrictspace=frozenseedingfock,
        localisometryselectionstrategy=frozenselectionbycount(3))

    globaldistillerspectrum = globaldistiller |> crystalspectrum
    # visualize(globaldistillerspectrum, title="Global Distiller")

    distillresult = distillation(globaldistillerspectrum, :courier => v -> abs(v) < 1e-5, :empty => v -> v > 1e-5, :filled => v -> v < -1e-5)

    courierseedingcenter::Offset = [2/3, 1/3] ∈ (blockedmodes |> getspace)
    courierseedingmodes::Subset{Mode} = circularregionmodes(courierseedingcenter, physicalmodes, 1.8)
    courierseedingregion::Subset{Offset} = Subset(m |> getpos for m in courierseedingmodes)
    # visualize(courierseedingregion, courierseedingcenter |> Subset, title="Courier Seeding Region", visualspace=euclidean(RealSpace, 2))
    courierseedingfock::FockSpace{Region} = FockSpace{Region}(courierseedingmodes)

    c3 = c6^2 |> recenter(courierseedingcenter)

    blockedcourierprojector = distillresult[:courier] |> crystalprojector
    blockedcouriercorrelation = idmap(blockedcourierprojector.outspace, blockedcourierprojector.outspace) - blockedcourierprojector

    localcourierseed = findlocalspstates(statecorrelations=blockedcouriercorrelation, regionfock=courierseedingfock, symmetry=c3, spectrumextractpredicate=v -> v < 5e-2)[1]
    fullcourierseed = localcourierseed + (c6 * localcourierseed.outspace) * localcourierseed * (c6 * localcourierseed.inspace)'

    crystalcourierseeds = crystalisometries(localisometry=fullcourierseed, crystalfock=blockedcorrelations.outspace, addinspacemomentuminfo=true)

    wanniercourierisometry = wannierprojection(
        crystalisometries=distillresult[:courier].eigenvectors, crystal=blockedcrystal, crystalseeds=crystalcourierseeds)

    # regionalrestriction(wanniercourierisometry, courierseedingfock) |> visualize

    couriercorrelations = wanniercourierisometry' * blockedcorrelations * wanniercourierisometry
    # couriercorrelations |> crystalspectrum |> visualize
    couriercorrelationspectrum = couriercorrelations |> crystalspectrum
    purifiedcorrelationspectrum = couriercorrelationspectrum |> roundingpurification
    # purifiedcorrelationspectrum |> visualize
    unpurecorrelations = couriercorrelations
    couriercorrelations = purifiedcorrelationspectrum |> FockMap

    blockedfilledprojector = distillresult[:filled] |> crystalprojector
    blockedfilledcorrelation = idmap(blockedfilledprojector.outspace, blockedfilledprojector.outspace) - blockedfilledprojector
    filledseed = findlocalspstates(statecorrelations=blockedfilledcorrelation, regionfock=frozenseedingfock, symmetry=c6, spectrumextractpredicate=v -> v < 1e-2, degeneracythreshold=1e-3)[3]

    crystalfilledseeds = crystalisometries(localisometry=filledseed, crystalfock=blockedcorrelations.outspace, addinspacemomentuminfo=true)
    wannierfilledisometry = wannierprojection(
        crystalisometries=distillresult[:filled].eigenvectors, crystal=blockedcrystal, crystalseeds=crystalfilledseeds)

    blockedemptyprojector = distillresult[:empty] |> crystalprojector
    blockedemptycorrelation = idmap(blockedemptyprojector.outspace, blockedemptyprojector.outspace) - blockedemptyprojector
    emptyseed = findlocalspstates(statecorrelations=blockedemptycorrelation, regionfock=frozenseedingfock, symmetry=c6, spectrumextractpredicate=v -> v < 1e-2, degeneracythreshold=1e-3)[3]

    crystalemptyseeds = crystalisometries(localisometry=emptyseed, crystalfock=blockedcorrelations.outspace, addinspacemomentuminfo=true)
    wannieremptyisometry = wannierprojection(
        crystalisometries=distillresult[:empty].eigenvectors, crystal=blockedcrystal, crystalseeds=crystalemptyseeds)

    filledcorrelations = wannierfilledisometry' * blockedcorrelations * wannierfilledisometry
    emptycorrelations = wannieremptyisometry' * blockedcorrelations * wannieremptyisometry

    isometry = wannierfilledisometry + wanniercourierisometry + wannieremptyisometry    # the total isometry (unitary)
    u_zipper = blocker' * isometry

    C_total = u_zipper' * correlations * u_zipper

    #modifiedTotalCorrelationMat = modified_filledcorrelations + correlations + modified_emptycorrelations
    
    swavefockspace = FockSpace(mode for mode in filledcorrelations |> getinspace if (mode|> getorbital(swave) |> rep) == [-1+0im])
    swavefockspace = FockSpace(swavefockspace, reflected=wannierfilledisometry |> getinspace |> getcrystal)
    swavefilledisometry = wannierfilledisometry[:, swavefockspace]
    
    swavefilledprojector = swavefilledisometry * swavefilledisometry'
    swavefilledcorrelations = idmap(swavefilledprojector |> getinspace) - swavefilledprojector
    filledcorrelations_tilde = wannierfilledisometry' * swavefilledcorrelations * wannierfilledisometry
    
    C_tilde = filledcorrelations_tilde + couriercorrelations + emptycorrelations 

    C_previous_tilde = u_zipper * C_tilde * u_zipper'
    C_previous= u_zipper * C_total * u_zipper'

    return Dict(
        :blocker => blocker,
        :blockedcorrelations => blockedcorrelations,
        :correlations => couriercorrelations, 
        :courierzipper => wanniercourierisometry' * blocker, 
        :filledzipper => wannierfilledisometry' * blocker,
        :emptyzipper => wannieremptyisometry' * blocker,
        :globaldistiller => globaldistiller, 
        :filledcorrelations => filledcorrelations, 
        :emptycorrelations => emptycorrelations,
        :entanglemententropy => entanglemententropy(couriercorrelationspectrum),
        :wanniercourierisometry => wanniercourierisometry,
        :wannierfilledisometry => wannierfilledisometry,
        :wannieremptyisometry => wannieremptyisometry, 
        :C_total => C_total,
        :filledcorrelations_tilde => filledcorrelations_tilde,
        :C_tilde => C_tilde,
        :isometry => isometry,
        :u_zipper => u_zipper,
        :C_previous_tilde => C_previous_tilde,
        :C_previous => C_previous,
        :frozenseedingfock => frozenseedingfock,
        :courierseedingfock => courierseedingfock,
        :distillresult => distillresult)
end 

# rg_tree = []
# function zer(correlations::FockMap)
#     i=0
#     push!(rg_tree, zipper(correlations))
#     try push!(rg_tree, zipper(correlations))
#         correlations = zipper(correlations)[:correlations]
#         i += 1
#     catch
# end

@time begin
rg1 = zipper(C)  
end  
@time begin
rg2 = zipper_even(rg1[:correlations]) 
end  
@time begin
rg3 = zipper(rg2[:correlations]) 
end  
@time begin
rg4 = zipper_even(rg3[:correlations])  
end

#####21/12/23###########################################################################################################################
# understand the current blocking function
function Base.:*(scale::Scale, crystalfock::FockSpace{Crystal})::FockMap
    crystal::Crystal = crystalfock |> getcrystal
    scaledcrystal::Crystal = scale * crystal
    unscaledblockedregion::Subset{Offset} = (scale |> inv) * scaledcrystal.unitcell
    bz::Subset{Momentum} = crystal |> brillouinzone
    basismodes::Subset{Mode} = crystalfock |> rep |> first
    scaledbz::Subset{Momentum} = scaledcrystal |> brillouinzone

    momentummappings::Base.Generator = (basispoint(scale * p) => p for p in bz)     #?
    mappingpartitions::Dict{Point, Vector{Point}} = foldl(momentummappings; init=Dict{Point, Vector{Point}}()) do d, (k, v)
        mergewith!(append!, d, LittleDict(k => [v]))
    end         #?

    ksubsets::Dict{Momentum, Subset{Mode}} = crystalfock |> crystalsubsets      # c(k) dictionary
    scaledfock::FockSpace = ((ksubsets[k] for k in mappingpartitions[kscaled]) |> subsetunion |> FockSpace
        for kscaled in scaledbz) |> fockspaceunion          # scaled Fock space

    # TODO: Use of latticeoff requires revisit, since latticeoff only truncates the decimals.
    unscaledblockedlatticeoffsets::Subset{Offset} = Subset(pbc(crystal, p) |> latticeoff for p in unscaledblockedregion)        # 
    unscaledblockedunitcellfock::FockSpace = spanoffset(basismodes, unscaledblockedlatticeoffsets)

    restrictedfourier::FockMap = fourier(bz, unscaledblockedunitcellfock)     # fourier(momentum points, FockSpace) gives a Fock Map representing a Fourier transform from the given FockSpace
    volumeratio::Number = (crystal |> vol) / (scaledcrystal |> vol)
    permutedfourier::FockMap = Zipper.permute(restrictedfourier, outspace=scaledfock) / sqrt(volumeratio)  # Rescale the FT due to space rescaling

    function repackfourierblocks(source::FockMap, kscaled::Momentum, partition::Subset{Mode})::FockMap
        partitionrows::FockMap = rows(source, partition |> FockSpace)
        inspace::FockSpace = (Subset(setattr(mode, :offset => kscaled, :pos => scale * convert(Point, mode))
            for mode in partitionrows.inspace |> orderedmodes)
            |> FockSpace)
        return FockMap(partitionrows.outspace, inspace, partitionrows |> rep)
    end

    repackedblocks::Base.Generator = (
        repackfourierblocks(permutedfourier, kscaled, partition)
        for (kscaled, partition) in Iterators.zip(scaledbz, scaledfock |> rep))
    blocking::FockMap = directsum(repackedblocks)
    return FockMap(blocking, inspace=FockSpace(blocking |> getinspace, reflected=scaledcrystal), outspace=crystalfock)'
end
#####19/12/23###########################################################################################################################
# testing for C

C_p = rg1[:filledzipper]' * rg1[:filledcorrelations] * rg1[:filledzipper] + rg1[:emptyzipper]' * rg1[:emptycorrelations] * rg1[:emptyzipper] + rg1[:courierzipper]' *  rg1[:correlations] * rg1[:courierzipper]  
Chern_number(groundstatespectrum(C_p |> crystalspectrum, perunitcellfillings=1))

# Just checking if doubling the bands gives anything weird. Answer: Really weird indeed!!!!! Revisit 'blocking'
C_doubled = blocked_correlation(C, 1)

C_doubled |> crystalspectrum |> visualize 

groundstate_doubled = groundstatespectrum(C_doubled |> crystalspectrum, perunitcellfillings=4) 
groundstate_doubled |> visualize

Chern_number(groundstate_doubled)


#####17/12/23###########################################################################################################################
C3_tilde = rg4[:filledzipper]' * rg4[:filledcorrelations_tilde] * rg4[:filledzipper] + rg4[:emptyzipper]' * rg4[:emptycorrelations] * rg4[:emptyzipper] + rg4[:courierzipper]' *  rg4[:correlations] * rg4[:courierzipper]
C2_tilde = rg3[:filledzipper]' * rg3[:filledcorrelations_tilde] * rg3[:filledzipper] + rg3[:emptyzipper]' * rg3[:emptycorrelations] * rg3[:emptyzipper] + rg3[:courierzipper]' *  C3_tilde * rg3[:courierzipper]
C1_tilde = rg2[:filledzipper]' * rg2[:filledcorrelations_tilde] * rg2[:filledzipper] + rg2[:emptyzipper]' * rg2[:emptycorrelations] * rg2[:emptyzipper] + rg2[:courierzipper]' *  C2_tilde * rg2[:courierzipper]

C_tilde = rg1[:filledzipper]' * rg1[:filledcorrelations_tilde] * rg1[:filledzipper] + rg1[:emptyzipper]' * rg1[:emptycorrelations] * rg1[:emptyzipper] + rg1[:courierzipper]' *  C1_tilde * rg1[:courierzipper]  

C_tilde |> crystalspectrum |> visualize

# N_bands = 2x4^d 

function blocked_correlation(correlations::FockMap, scaling::Int64)
    crystalfockspace = correlations |> getoutspace
    scale = Scale([2 0; 0 2]^scaling, crystalfockspace |> getcrystal |> getspace)
    blockresult = blocking(:action => scale, :correlations => correlations, :crystal => crystalfockspace |> getcrystal)
    return blockresult[:correlations]
end

C_tilde_sym = blocked_correlation(C_tilde, 4)

C_tilde_sym |> crystalspectrum |> visualize

groundstate_sym = groundstatespectrum(C_tilde_sym |> crystalspectrum, perunitcellfillings=86) 
groundstate_sym|> visualize

Chern_number(groundstate_sym)

# Just checking if doubling the bands gives anything weird. Answer: No. Wait. Really weird indeed!!!!!
C_p = rg1[:filledzipper]' * rg1[:filledcorrelations] * rg1[:filledzipper] + rg1[:emptyzipper]' * rg1[:emptycorrelations] * rg1[:emptyzipper] + rg1[:courierzipper]' *  rg1[:correlations] * rg1[:courierzipper]  
Chern_number(groundstatespectrum(C_p |> crystalspectrum, perunitcellfillings=1))

C_doubled = blocked_correlation(C_p, 4)

C_doubled |> crystalspectrum |> visualize 

groundstate_doubled = groundstatespectrum(C_doubled |> crystalspectrum, perunitcellfillings=257) 
groundstate_doubled |> visualize

Chern_number(groundstate_doubled)

Chern_number(groundstatespectrum(C |> crystalspectrum, perunitcellfillings=1))

#####08/12/23###########################################################################################################################
# get real space representation of correlation matrix
# fourier_trans = fourier(C |> getinspace, sparsefock(modes, crystal |> latticepoints))

# C_real = fourier_trans' * C * fourier_trans


# mA = Mode(:offset=>(0,0)∈triangular, :pos=>(1/3,2/3)∈ triangular, :flavor=>1)
# mB = Mode(:offset=>(0,0)∈triangular, :pos=>(2/3,1/3)∈ triangular, :flavor=>1)

# C_real[mA, mB]
# bonds

# Construct an effective Hamiltonian that can be more physical (less flat band) 


#####05/12/23###########################################################################################################################
C3_tilde = rg4[:filledzipper]' * rg4[:filledcorrelations_tilde] * rg4[:filledzipper] + rg4[:emptyzipper]' * rg4[:emptycorrelations] * rg4[:emptyzipper] + rg4[:courierzipper]' *  rg4[:correlations] * rg4[:courierzipper]
C2_tilde = rg3[:filledzipper]' * rg3[:filledcorrelations_tilde] * rg3[:filledzipper] + rg3[:emptyzipper]' * rg3[:emptycorrelations] * rg3[:emptyzipper] + rg3[:courierzipper]' *  C3_tilde * rg3[:courierzipper]
C1_tilde = rg2[:filledzipper]' * rg2[:filledcorrelations_tilde] * rg2[:filledzipper] + rg2[:emptyzipper]' * rg2[:emptycorrelations] * rg2[:emptyzipper] + rg2[:courierzipper]' *  C2_tilde * rg2[:courierzipper]

C_tilde = rg1[:filledzipper]' * rg1[:filledcorrelations_tilde] * rg1[:filledzipper] + rg1[:emptyzipper]' * rg1[:emptycorrelations] * rg1[:emptyzipper] + rg1[:courierzipper]' *  C1_tilde * rg1[:courierzipper]  

C_tilde |> crystalspectrum |>  visualize

C_tilde_mat = C_tilde |> rep
col_list = []
for col in range(1, size(C_tilde_mat[1,:])[1])
    if C_tilde_mat[1,col] != 0
        push!(col_list, col)
    end
end
col_list
size(col_list)
# N_bands = 2x4^d 

function blocking_scale(correlations::FockMap, scaling::Int64)
    crystalfockspace = correlations |> getoutspace
    scale = Scale([2 0; 0 2]^scaling, crystalfockspace |> getcrystal |> getspace)
    blockresult = blocking(:action => scale, :correlations => correlations, :crystal => crystalfockspace |> getcrystal)
    return blockresult
end

blockresult = blocking_scale(C_tilde, 4)
# blockedcrystal::Crystal = blockresult[:crystal]
C_tilde_sym = blockresult[:correlations]
# blocker = blockresult[:transformer]

C_tilde_sym_mat = C_tilde_sym |> rep
col_sym_list = []
for col in range(1, size(C_tilde_sym_mat[1,:])[1])
    if C_tilde_sym_mat[1,col] != 0
        push!(col_sym_list, col)
    end
end
col_sym_list
size(col_sym_list)

C_tilde_sym |> crystalspectrum |> visualize

# groundstatespectrum(C |> crystalspectrum, perunitcellfillings=2) |> visualize

groundstate_sym = groundstatespectrum(C_tilde_sym |> crystalspectrum, perunitcellfillings=86) 
groundstate_sym|> visualize

get_Chern_number(groundstate_sym)


get_Chern_number(groundstatespectrum(C |> crystalspectrum, perunitcellfillings=1))

# H_eff = correlation
H_eff_tilde = groundstate_sym |> FockMap
# H_eff |> crystalspectrum |> visualize
# blocker' * H_eff_tilde * blocker |> crystalspectrum |> visualize


# blocker' * H_eff_tilde * blocker |> crystalspectrum |> visualize
using Plots
Clist = [0.997,0.496,0.169,0.039, 0.010]
dlist = [0,1,2,3,4]
plot!(dlist, Clist, xlabel = 'd', ylabel = 'C', legend = false)

#####04/12/23###########################################################################################################################
C3 = rg3[:correlations]
C2 = rg3[:filledzipper]' * rg3[:filledcorrelations] * rg3[:filledzipper] + rg3[:emptyzipper]' * rg3[:emptycorrelations] * rg3[:emptyzipper] + rg3[:courierzipper]' *  rg3[:correlations] * rg3[:courierzipper]
C2 |> visualize
C1 = rg2[:filledzipper]' * rg2[:filledcorrelations] * rg2[:filledzipper] + rg2[:emptyzipper]' * rg2[:emptycorrelations] * rg2[:emptyzipper] + rg2[:courierzipper]' *  C2 * rg2[:courierzipper]
C1 |> crystalspectrum |> visualize
C0 = rg1[:filledzipper]' * rg1[:filledcorrelations] * rg1[:filledzipper] + rg1[:emptyzipper]' * rg1[:emptycorrelations] * rg1[:emptyzipper] + rg1[:courierzipper]' *  rg1[:correlations] * rg1[:courierzipper]  
get_Chern_number(groundstatespectrum(rg3[:correlations]|>crystalspectrum, perunitcellfillings=1))

C3_tilde = rg3[:correlations]#rg4[:filledzipper]' * rg4[:filledcorrelations_tilde] * rg4[:filledzipper] + rg4[:emptyzipper]' * rg4[:emptycorrelations] * rg4[:emptyzipper] + rg4[:courierzipper]' *  rg4[:correlations] * rg4[:courierzipper]
C2_tilde = rg3[:filledzipper]' * rg3[:filledcorrelations_tilde] * rg3[:filledzipper] + rg3[:emptyzipper]' * rg3[:emptycorrelations] * rg3[:emptyzipper] + rg3[:courierzipper]' *  C3_tilde * rg3[:courierzipper]
C1_tilde = rg2[:filledzipper]' * rg2[:filledcorrelations_tilde] * rg2[:filledzipper] + rg2[:emptyzipper]' * rg2[:emptycorrelations] * rg2[:emptyzipper] + rg2[:courierzipper]' *  C2_tilde * rg2[:courierzipper]
C_tilde = rg1[:filledzipper]' * rg1[:filledcorrelations_tilde] * rg1[:filledzipper] + rg1[:emptyzipper]' * rg1[:emptycorrelations] * rg1[:emptyzipper] + rg1[:courierzipper]' *  C1_tilde * rg1[:courierzipper]  

C_tilde |> crystalspectrum |>  visualize
groundstateprime = groundstatespectrum(C_tilde|> crystalspectrum, perunitcellfillings=1.0)  |> visualize

C_tilde

rg1[:courierzipper]' *  rg1[:correlations] * rg1[:courierzipper] |> visualize

groundstatespectrum(rg1[:courierzipper]' *  rg1[:correlations] * rg1[:courierzipper] |>crystalspectrum, perunitcellfillings=1) |> visualize
get_Chern_number(groundstatespectrum(rg1[:courierzipper]' *  rg1[:correlations] * rg1[:courierzipper] |>crystalspectrum, perunitcellfillings=1))


C_tilde_mat = C_tilde |> rep
col_list = []
for col in range(1, size(C_tilde_mat[1,:])[1])
    if C_tilde_mat[1,col] != 0
        push!(col_list, col)
    end
end
col_list
size(col_list)
128/2==4^3
[2 0; 0 2]^3

# function get_blocked_correlation(correlations::FockMap)::
    
# end

crystalfockspace = C_tilde |> getoutspace
scale = Scale([2 0; 0 2]^3, crystalfockspace |> getcrystal |> getspace)
blockresult = blocking(:action => scale, :correlations => C_tilde, :crystal => crystalfockspace |> getcrystal)

blockedcrystal::Crystal = blockresult[:crystal]
correct_correlations = blockresult[:correlations]
# blocker = blockresult[:transformer]

correct_groundstates = groundstatespectrum(correct_correlations |> crystalspectrum, perunitcellfillings=1) |> visualize

get_Chern_number(correct_groundstates)

# newcrystalspectrum = crystalspectrum(correction |> crystalsubmaps, crystal=Crystal(unitcell, [12, 12]))

newcrystalspectrum |> visualize

rg3[:correlations] |> getoutspace |> getcrystal

Cmat = C|>rep
# Effective Hamiltonian
# H_eff = -log(C.inversed - 1)

# C |> getinverse
# effectiveHamiltonian = -log(inv(C)-idmap(C.outspace))

function effectiveHamiltonianSpectrum(momentumfockmaps; crystal::Crystal)::CrystalSpectrum
    crystaleigenmodes::Dict{Momentum, Subset{Mode}} = Dict()
    crystaleigenvalues::Dict{Mode, Number} = Dict()
    H_eigenvalues::Dict{Mode, Number} = Dict()
    crystaleigenvectors::Dict{Momentum, FockMap} = Dict()
    for (k, fockmap) in momentumfockmaps
        eigenspectrum::EigenSpectrum = eigspech(fockmap, :offset => k)
        crystaleigenmodes[k] = Subset(m for (m, _) in eigenspectrum |> geteigenvalues)
        crystaleigenvectors[k] = eigenspectrum |> geteigenvectors
        for (m, v) in eigenspectrum |> geteigenvalues
            crystaleigenvalues[m] = v
            H_eigenvalues[k] =  -log(1/v-1)
        end
    end
    return CrystalSpectrum(crystal, crystaleigenmodes, H_eigenvalues, crystaleigenvectors)
end

effectiveHamiltonianSpectrum(fockmap::FockMap)::CrystalSpectrum = crystalspectrum(fockmap |> crystalsubmaps, crystal=fockmap.inspace |> getcrystal)

groundstatespectrumnew = groundstatespectrum(C |> effectiveHamiltonianSpectrum, perunitcellfillings=1.0) 
groundstatespectrumnew |> visualize

###
C_tilde |> crystalspectrum |> visualize
groundstatemm = groundstatespectrum(correct_correlations |> crystalspectrum, perunitcellfillings=1.0/3) 
groundstatemm |> visualize

##########################################################################################
correct_C_tilde_mat = correct_correlations |> rep  #
visualize(correct_correlations, colrange=1:1000, rowrange=1:1000)

correct_C_tilde_mat[1,2305]

colu_list = []
for col in range(1, size(correct_C_tilde_mat[1,:])[1])
    if C_tilde_mat[1,col] != 0
        push!(colu_list, col)
    end
end
size(colu_list)



#####30/11/23###########################################################################################################################
correlations = C
crystalfock = C.outspace

scale = Scale([2 0; 0 2], crystalfock |> getcrystal |> getspace)
blockresult = blocking(:action => scale, :correlations => correlations, :crystal => crystalfock |> getcrystal)

blockedcrystal::Crystal = blockresult[:crystal]
blockedcorrelations::FockMap = blockresult[:correlations]
blocker = blockresult[:transformer]

function circularregionmodes(origin::Offset, physicalmodes::Subset{Mode}, radius::Number)::Subset{Mode}
    currentspace::RealSpace = correlations.outspace |> getcrystal |> getspace |> orthogonalspace
    physicalnorm = m -> lineartransform(currentspace, m |> getpos) |> norm
    return filter(p -> physicalnorm(p - origin) < radius, physicalmodes)
end

crystalpoints::Subset{Offset} = latticepoints(blockedcrystal)
samplepoints::Subset{Offset} = crystalpoints + c6^2 * crystalpoints + c6^4 * crystalpoints
blockedmodes::Subset{Mode} = quantize(:pos, blockedcrystal.unitcell, 1)
physicalmodes::Subset{Mode} = spanoffset(blockedmodes, samplepoints)


frozenseedingmodes::Subset{Mode} = circularregionmodes(triangular |> getorigin, physicalmodes, 2)
# visualize(frozenseedingregion, title="Frozen Seeding Region", visualspace=euclidean(RealSpace, 2))
frozenseedingfock::FockSpace = FockSpace{Region}(frozenseedingmodes)

regioncorrelations(blockedcorrelations, frozenseedingfock) |> eigspech |> visualize
plot(scatter(y=matr,ms=2, ma=0.5)) 

globaldistiller = globaldistillerhamiltonian(
    correlations=blockresult[:correlations],
    restrictspace=frozenseedingfock,
    localisometryselectionstrategy=frozenselectionbycount(3))

globaldistillerspectrum = globaldistiller |> crystalspectrum
# visualize(globaldistillerspectrum, title="Global Distiller")

distillresult = distillation(globaldistillerspectrum, :courier => v -> abs(v) < 1e-5, :empty => v -> v > 1e-5, :filled => v -> v < -1e-5)


#####29/11/23###########################################################################################################################
scale = Scale([2 0; 0 2], triangular)
blocker = scale * C.outspace
H = energyspectrum |> FockMap
blockedH = blocker * H * blocker'
filledP = rg1[:wannierfilledisometry] * rg1[:wannierfilledisometry]'
physfilledH = blocker' * filledP * blockedH * filledP * blocker
physfilledH |> crystalspectrum |> visualize

filledP = rg1[:wanniercourierisometry] * rg1[:wanniercourierisometry]' 
physfilledH =  filledP * blockedH * filledP 
physfilledH |> crystalspectrum |> visualize

rg1[:wannieremptyisometry] * rg1[:wannieremptyisometry]'

#####27/11/23###########################################################################################################################
C3_tilde = rg4[:filledzipper]' * rg4[:filledcorrelations_tilde] * rg4[:filledzipper] + rg4[:emptyzipper]' * rg4[:emptycorrelations] * rg4[:emptyzipper] + rg4[:courierzipper]' *  rg4[:correlations] * rg4[:courierzipper]
C2_tilde = rg3[:filledzipper]' * rg3[:filledcorrelations_tilde] * rg3[:filledzipper] + rg3[:emptyzipper]' * rg3[:emptycorrelations] * rg3[:emptyzipper] + rg3[:courierzipper]' *  C3_tilde * rg3[:courierzipper]
C1_tilde = rg2[:filledzipper]' * rg2[:filledcorrelations_tilde] * rg2[:filledzipper] + rg2[:emptyzipper]' * rg2[:emptycorrelations] * rg2[:emptyzipper] + rg2[:courierzipper]' *  C2_tilde * rg2[:courierzipper]
C_tilde = rg1[:filledzipper]' * rg1[:filledcorrelations_tilde] * rg1[:filledzipper] + rg1[:emptyzipper]' * rg1[:emptycorrelations] * rg1[:emptyzipper] + rg1[:courierzipper]' *  C1_tilde * rg1[:courierzipper]  

# Create effective Hamiltonian
# function effectiveHamiltonian(correlations::FockMap, effective_crystal::Crystal)
#     jjj
#     ddd
# end

#####26/11/23###########################################################################################################################
C3 = rg3[:correlations]
C2 = rg3[:filledzipper]' * rg3[:filledcorrelations] * rg3[:filledzipper] + rg3[:emptyzipper]' * rg3[:emptycorrelations] * rg3[:emptyzipper] + rg3[:courierzipper]' *  rg3[:correlations] * rg3[:courierzipper]
C2 |> visualize
C1 = rg2[:filledzipper]' * rg2[:filledcorrelations] * rg2[:filledzipper] + rg2[:emptyzipper]' * rg2[:emptycorrelations] * rg2[:emptyzipper] + rg2[:courierzipper]' *  C2 * rg2[:courierzipper]
C1 |> crystalspectrum |> visualize
C0 = rg1[:filledzipper]' * rg1[:filledcorrelations] * rg1[:filledzipper] + rg1[:emptyzipper]' * rg1[:emptycorrelations] * rg1[:emptyzipper] + rg1[:courierzipper]' *  rg1[:correlations] * rg1[:courierzipper]  
visualize(C0, colrange=1:1000, rowrange=1:1000) 
C0 |> crystalspectrum |> visualize
get_Chern_number(C0 |> crystalspectrum)
get_Chern_number(rg1[:courierzipper]' *  rg1[:correlations] * rg1[:courierzipper] |>crystalspectrum )
get_Chern_number(C|> crystalspectrum)

C3_tilde = rg3[:correlations] #rg4[:filledzipper]' * rg4[:filledcorrelations_tilde] * rg4[:filledzipper] + rg4[:emptyzipper]' * rg4[:emptycorrelations] * rg4[:emptyzipper] + rg4[:courierzipper]' *  rg4[:correlations] * rg4[:courierzipper]
C2_tilde = rg3[:filledzipper]' * rg3[:filledcorrelations_tilde] * rg3[:filledzipper] + rg3[:emptyzipper]' * rg3[:emptycorrelations] * rg3[:emptyzipper] + rg3[:courierzipper]' *  C3_tilde * rg3[:courierzipper]
C1_tilde = rg2[:filledzipper]' * rg2[:filledcorrelations_tilde] * rg2[:filledzipper] + rg2[:emptyzipper]' * rg2[:emptycorrelations] * rg2[:emptyzipper] + rg2[:courierzipper]' *  C2_tilde * rg2[:courierzipper]
C_tilde = rg1[:filledzipper]' * rg1[:filledcorrelations_tilde] * rg1[:filledzipper] + rg1[:emptyzipper]' * rg1[:emptycorrelations] * rg1[:emptyzipper] + rg1[:courierzipper]' *  C1_tilde * rg1[:courierzipper]  

state_tilde = C_tilde |> crystalspectrum
state_tilde |> geteigenvectors
get_Chern_number(state_tilde)
get_Chern_number(groundstates)

# state_tilde2 = C2_tilde |> crystalspectrum
# get_Chern_number(state_tilde2) need to fix for general lattice

#####25/11/23###########################################################################################################################
test_corr = rg4[:filledzipper]' * rg4[:filledcorrelations] * rg4[:filledzipper] + rg4[:courierzipper]' * rg4[:correlations] * rg4[:courierzipper] + rg4[:emptyzipper]' * rg4[:emptycorrelations] * rg4[:emptyzipper]
rg4[:C_previous] - test_corr |> visualize # valid or not???
# anyway, try it
C3_tilde = rg4[:filledzipper]' * rg4[:filledcorrelations_tilde] * rg4[:filledzipper] + rg4[:emptyzipper]' * rg4[:emptycorrelations] * rg4[:emptyzipper] + rg4[:courierzipper]' *  rg4[:correlations] * rg4[:courierzipper]
C2_tilde = rg3[:filledzipper]' * rg3[:filledcorrelations_tilde] * rg3[:filledzipper] + rg3[:emptyzipper]' * rg3[:emptycorrelations] * rg3[:emptyzipper] + rg3[:courierzipper]' *  C3_tilde * rg3[:courierzipper]
C1_tilde = rg2[:filledzipper]' * rg2[:filledcorrelations_tilde] * rg2[:filledzipper] + rg2[:emptyzipper]' * rg2[:emptycorrelations] * rg2[:emptyzipper] + rg2[:courierzipper]' *  C2_tilde * rg2[:courierzipper]
C_tilde = rg1[:filledzipper]' * rg1[:filledcorrelations_tilde] * rg1[:filledzipper] + rg1[:emptyzipper]' * rg1[:emptycorrelations] * rg1[:emptyzipper] + rg1[:courierzipper]' *  C1_tilde * rg1[:courierzipper]  
C_tilde_mat = C_tilde |> rep 
size(C_tilde_mat[:,1])[1] 
C_tilde_mat[1,3] == 0
n = 0   # count numbers of momenta couplings
col_list = []
for col in range(1, size(C_tilde_mat[1,:])[1])
    if C_tilde_mat[1,col] != 0
        push!(col_list, col)
    end
end
col_list
size(col_list)
512==2^9 

# how to change color to plot the whole BZ
using Plots; gr()
using ColorSchemes, Colors
function visualizeFM(fockmap::FockMap; title::String = "", rowrange = :, colrange = :)
    source = rep(fockmap)[rowrange, colrange]
    layout = Layout(yaxis=attr(autorange="reversed"))
    subplots = [plot(heatmap(z=real(source)), layout) plot(heatmap(z=imag(source)), layout)]
    relayout!(subplots, title_text=title)
    subplots
end

visualizeFM(C_tilde, colrange=1:1000, rowrange=1:1000)
# , cmax=1.0, cmin= -1.0, , zmin=-1.0, zmax=1.0

# Chern number
function computequantumgeometrictensor2d(state::CrystalSpectrum)
    statevectors::Dict{Momentum, FockMap} = state|>geteigenvectors
    Δkx = (1 / size(crystal)[1], 0) ∈ kspace
    Δky = (0, 1 / size(crystal)[2]) ∈ kspace

    kprojectors::Dict{Momentum, FockMap} = Dict(k => kstate * kstate' for (k, kstate) in statevectors)
    kcorrelations::Dict{Momentum, FockMap} = Dict(k => idmap(projector|>getoutspace) - projector for (k, projector) in kprojectors)
    ΔCkxs::Dict{Momentum, FockMap} = Dict(k => inplaceadd(-kcorr, kcorrelations[k + Δkx |> basispoint]) for (k, kcorr) in kcorrelations)
    ΔCkys::Dict{Momentum, FockMap} = Dict(k => inplaceadd(-kcorr, kcorrelations[k + Δky |> basispoint]) for (k, kcorr) in kcorrelations)
    
    function kqgtensor(k::Momentum)::FockMap
        gxx = statevectors[k]' * (ΔCkxs[k]' * kcorrelations[k] * ΔCkxs[k]) * statevectors[k]
        gxy = statevectors[k]' * (ΔCkxs[k]' * kcorrelations[k] * ΔCkys[k]) * statevectors[k]
        gyx = statevectors[k]' * (ΔCkys[k]' * kcorrelations[k] * ΔCkxs[k]) * statevectors[k]
        gyy = statevectors[k]' * (ΔCkys[k]' * kcorrelations[k] * ΔCkys[k]) * statevectors[k]

        xmode::Mode = Mode(:offset=>k, :axis=>(:x))
        ymode::Mode = Mode(:offset=>k, :axis=>(:y))

        fockspace::FockSpace = FockSpace([xmode, ymode])
        return FockMap(fockspace, fockspace, [rep(gxx)[1,1] rep(gxy)[1,1]; rep(gyx)[1,1] rep(gyy)[1,1]])
    end

    return Dict(k=>k|>kqgtensor for (k, _) in statevectors)
end

Base.:real(fockmap::FockMap)::FockMap = FockMap(fockmap|>getoutspace, fockmap|>getinspace, fockmap|>rep|>real)
Base.:imag(fockmap::FockMap)::FockMap = FockMap(fockmap|>getoutspace, fockmap|>getinspace, fockmap|>rep|>imag)

kqgtensors = computequantumgeometrictensor2d(groundstates)
gxx = directsum(tensor[1,2] for (k, tensor) in kqgtensors)
barecrystal = Crystal(triangular|>getorigin, crystal|>size)
barecrystalfock = FockSpace(gxx|>getoutspace, reflected=barecrystal)
gxx = FockMap(gxx, inspace=barecrystalfock, outspace=barecrystalfock, performpermute=false)|>real
gxx|>crystalspectrum|>visualize

bc = directsum(tensor[1,2]|>imag for (k, tensor) in kqgtensors)
barecrystalfock = FockSpace(bc|>getoutspace, reflected=barecrystal)
berrycurvature = 2 * FockMap(bc, inspace=barecrystalfock, outspace=barecrystalfock, performpermute=false)

(berrycurvature|>rep|>sum) / (2π)

berrycurvature |>visualize

#####23/11/23###########################################################################################################################
# function UV_Correlations(IR_correlation::FockMap)::FockMap
#     fff
# end


# wrong, need to add filled and empty terms
C3_tilde = rg4[:C_previous_tilde]
C2_tilde = rg3[:courierzipper]' *  C3_tilde * rg3[:courierzipper]
C1_tilde = rg2[:courierzipper]' *  C2_tilde * rg2[:courierzipper]
C_tilde = rg1[:courierzipper]' *  C1_tilde * rg1[:courierzipper]  
C
# mat = C_tilde |> rep |> Matrix |> Hermitian  
# plot(scatter(y=mat |> eigvals)) 

C_tilde |> visualize 
visualize(C_tilde, colrange=1:1000, rowrange=1:1000)

#####22/11/23###########################################################################################################################
"""
rg1filled = rg1[:filledcorrelations] |> rep
rg2filled = rg2[:filledcorrelations] |> rep

filled_Fock1 = rg1[:filledcorrelations] |> getinspace 
list1 = [] 
for mode in filled_Fock1
    push!(list1, mode |> getorbital(swave))
    # getorbital(swave)
    # if mode|> getorbital(swave) == swave
    #     push!(list, mode)
    # end
end 

filled_Fock2 = rg2[:filledcorrelations] |> getinspace 
list2 = [] 
for mode in filled_Fock2
    push!(list2, mode |> getorbital(swave))
    #push!(list2, mode |> getorbital(swave))
    # getorbital(swave)
    # if mode|> getorbital(swave) == swave
    #     push!(list, mode)
    # end
end 

list2p = [] 
for mode in filled_Fock2
    if mode|> getorbital(swave) |> rep == [-1+0im]
        push!(list2p, mode)
    end
    #push!(list2, mode |> getorbital(swave))
    # getorbital(swave)
    # if mode|> getorbital(swave) == swave
    #     push!(list, mode)
    # end
end 

swavefockspace = FockSpace(mode for mode in rg2[:filledcorrelations] |> getinspace if (mode|> getorbital(swave) |> rep) == [-1+0im])
swavefockspace = FockSpace(swavefockspace, reflected=rg2[:wannierfilledisometry] |> getinspace |> getcrystal)
swavefilledisometry = rg2[:wannierfilledisometry][:, swavefockspace]

swavefilledprojector = swavefilledisometry * swavefilledisometry'
swavefilledcorrelations = idmap(swavefilledprojector |> getinspace) - swavefilledprojector
filledcorrelations_tilde = rg2[:wannierfilledisometry]' * swavefilledcorrelations * rg2[:wannierfilledisometry]
filledcorrelations_tilde |> visualize

C_tilde = filledcorrelations_tilde + rg2[:correlations] + rg2[:emptycorrelations] 
C_tilde |> visualize
C_previous_tilde = rg2[:u_zipper] * C_tilde * rg2[:u_zipper]'  
C_previous_tilde |> visualize
mat = C_previous_tilde |> rep |> Matrix |> Hermitian
plot(scatter(y=mat |> eigvals))

filled_Fock3 = rg3[:filledcorrelations] |> getinspace 
list3 = [] 
for mode in filled_Fock3
    if mode|> getorbital(swave) == swaveminus
        push!(list, mode)
    end
    #push!(list3, mode |> getorbital(swave))
    # getorbital(swave)
    # if mode|> getorbital(swave) == swave
    #     push!(list, mode)
    # end
end 

filled_Fock4 = rg4[:filledcorrelations] |> getinspace 
list4 = [] 
for mode in filled_Fock4
    push!(list4, mode |> getorbital(swave))
    # getorbital(swave)
    # if mode|> getorbital(swave) == swave
    #     push!(list, mode)
    # end
end 

mode for mode in filledcorrelations |> getinspace if swave == (mode |> getorbital(swave))

rg1[:C_previous] |> visualize

crystal |> sitepoints |> visualize

rg3[:correlations] |> getoutspace |> getcrystal |> sitepoints |> visualize

# rg1[:filledzipper]' * rg1[:filledcorrelations_tilde] * rg1[:filledzipper] |> visualize
# rg1[:emptyzipper]' * rg1[:emptycorrelations] * rg1[:emptyzipper] |> visualize
# rg1[:filledzipper]' * rg1[:filledcorrelations] * rg1[:filledzipper]+ rg1[:emptyzipper]' * rg1[:emptycorrelations] * rg1[:emptyzipper] + rg1[:courierzipper]' * rg1[:correlations] * rg1[:courierzipper] |> visualize
# rg1[:courierzipper]' * rg1[:correlations] * rg1[:courierzipper] |> visualize

# corr = rg1[:filledcorrelations_tilde] + rg1[:correlations] + rg1[:emptycorrelations]
# rg1[:u_zipper]
# physcorr |> visualize
# physcorr = rg1[:u_zipper] * corr * rg1[:u_zipper]'
# PhyscorrMat = physcorr |> rep |> Matrix
# import Pkg
# using Plots
# layout = Layout(yaxis=attr(autorange="reversed"))
# plot(heatmap(z=real(PhyscorrMat)), layout)

plot(scatter(y=physcorr |> rep |> Matrix |> Hermitian |> eigvals))

rg1[:C_previous_tilde] |> visualize

mat = rg1[:C_previous_tilde] |> rep |> Matrix |> Hermitian

###########################################################################################################################

#####22/11/23#####################################################################################
crystalfock = C.outspace     # map to the physical space

scale = Scale([2 0; 0 2], crystalfock |> getcrystal |> getspace)    
blockresult = blocking(:action => scale, :correlations => C, :crystal => crystalfock |> getcrystal)      
    
blockedcrystal::Crystal = blockresult[:crystal]

crystalpoints::Subset{Offset} = latticepoints(blockedcrystal)
blockedmodes::Subset{Mode} = quantize(:pos, blockedcrystal.unitcell, 1)
physicalmodes::Subset{Mode} = spanoffset(blockedmodes, crystalpoints)

function circularregionmodes(origin::Offset, physicalmodes::Subset{Mode}, radius::Number)::Subset{Mode}
    currentspace::RealSpace = C.outspace |> getcrystal |> getspace |> orthogonalspace
    physicalnorm = m -> lineartransform(currentspace, m |> getpos) |> norm
    return filter(p -> physicalnorm(p - origin) < radius, physicalmodes)
end

frozenseedingmodes::Subset{Mode} = circularregionmodes(triangular |> getorigin, physicalmodes, 4)         # specify modes in the distill region
    # frozenseedingregion::Subset{Offset} = Subset(m |> getpos for m in frozenseedingmodes)
    # visualize(frozenseedingregion, title="Frozen Seeding Region", visualspace=euclidean(RealSpace, 2))
frozenseedingfock::FockSpace = FockSpace{Region}(frozenseedingmodes)

# filled2physical = rg1[:courierzipper]' * rg2[:filledzipper]' * rg2[:filledcorrelations] * rg2[:filledzipper] * rg1[:courierzipper]

function visualize(state::RegionState{2}; title::String = "", markersizemultiplier::Real = 120)
    function generatestateplot(spstate::FockMap)
        spmode::Mode = spstate |> getinspace |> first
        columnspectrum::Base.Generator = (m => spstate[m, spmode] for m in spstate |> getoutspace |> orderedmodes)
        positions::Vector{Offset} = [v.first |> getpos for v in columnspectrum]
        mesh::Matrix{Float64} = hcat(map(p -> p |> euclidean |> vec, positions)...)
        markersizes::Vector{Float64} = [v.second |> abs for v in columnspectrum]
        normalizedmarkersizes::Vector{Float64} = markersizes / norm(markersizes) * markersizemultiplier
        markercolors::Vector = [convert(RGB{Float64}, HSV(angle(v.second) / 2π * 360, 1, 1)) for v in columnspectrum]
        return scatter(
            x=mesh[1, :], y=mesh[2, :], mode="markers",
            marker=attr(
                symbol="circle",
                size=normalizedmarkersizes,
                color=markercolors))
    end

    scatters::Vector = [spstate |> generatestateplot for (_, spstate) in state]
    subplotnames::Base.Generator = (mode |> string for mode in state |> getinspace |> orderedmodes |> indexmodes)
    fig = make_subplots(rows=1, cols=scatters |> length, subplot_titles=hcat(subplotnames...))
    for (n, scatter) in enumerate(scatters)
        add_trace!(fig, scatter, row=1, col=n)
    end
    relayout!(fig, title_text=title)
    fig
end

regionalrestriction(rg1[:wannierfilledisometry], frozenseedingfock) |> visualize

rg2filled_physical_isometry = rg1[:wanniercourierisometry]*rg2[:blocker]'*rg2[:wannierfilledisometry]
regionalrestriction(rg2filled_physical_isometry, frozenseedingfock) |> visualize

rg3filled_physical_isometry = rg1[:wanniercourierisometry]*rg2[:blocker]'*rg2[:wanniercourierisometry]*rg3[:blocker]'*rg3[:wannierfilledisometry]
regionalrestriction(rg3filled_physical_isometry, frozenseedingfock) |> visualize

##################################################################################################################################



c6.transformmatrix |> eigvals

pplus = findeigenfunction(c6, eigenvalue=(exp(im * π/3)))

swavefilledfockspace = FockSpace(mode for mode in wannierfilledunitary |> getinspace if (mode |> getorbital(swave)) == swave)

correlations[swavefilledfockspace, swavefilledfockspace]




H = energyspectrum |> FockMap
energyspectrum |> visualize

rg1[:blocker] * H * rg1[:blocker]' |> crystalspectrum |> visualize

rg1[:courierzipper] * H * rg1[:courierzipper]' |> crystalspectrum |> visualize
rg1[:filledzipper] * H * rg1[:filledzipper]' |> crystalspectrum |> visualize
rg1[:emptyzipper] * H * rg1[:emptyzipper]' |> crystalspectrum |> visualize

courierH = rg1[:courierzipper] * H * rg1[:courierzipper]'
rg2[:blocker] * courierH * rg2[:blocker]' |> crystalspectrum |> visualize

rg2[:courierzipper] * rg1[:courierzipper] * H * rg1[:courierzipper]' * rg2[:courierzipper]' |> crystalspectrum |> visualize
rg2[:filledzipper] * rg1[:courierzipper] * H * rg1[:courierzipper]' * rg2[:filledzipper]' |> crystalspectrum |> visualize
rg2[:emptyzipper] * rg1[:courierzipper] * H * rg1[:courierzipper]' * rg2[:emptyzipper]' |> crystalspectrum |> visualize

courierH = rg2[:courierzipper] * courierH * rg2[:courierzipper]'
rg3[:blocker] * courierH * rg3[:blocker]' |> crystalspectrum |> visualize

rg3[:courierzipper] * rg2[:courierzipper] * rg1[:courierzipper] * H * rg1[:courierzipper]' * rg2[:courierzipper]' * rg3[:courierzipper]' |> crystalspectrum |> visualize
rg3[:filledzipper] * rg2[:courierzipper] * rg1[:courierzipper] * H * rg1[:courierzipper]' * rg2[:courierzipper]' * rg3[:filledzipper]' |> crystalspectrum |> visualize
rg3[:emptyzipper] * rg2[:courierzipper] * rg1[:courierzipper] * H * rg1[:courierzipper]' * rg2[:courierzipper]' * rg3[:emptyzipper]' |> crystalspectrum |> visualize

courierH = rg3[:courierzipper] * courierH * rg3[:courierzipper]'
rg4[:blocker] * courierH * rg4[:blocker]' |> crystalspectrum |> visualize

rg4[:courierzipper] * rg3[:courierzipper] * rg2[:courierzipper] * rg1[:courierzipper] * H * rg1[:courierzipper]' * rg2[:courierzipper]' * rg3[:courierzipper]' * rg4[:courierzipper]' |> crystalspectrum |> visualize
rg4[:filledzipper] * rg3[:courierzipper] * rg2[:courierzipper] * rg1[:courierzipper] * H * rg1[:courierzipper]' * rg2[:courierzipper]' * rg3[:courierzipper]' * rg4[:filledzipper]' |> crystalspectrum |> visualize
rg4[:emptyzipper] * rg3[:courierzipper] * rg2[:courierzipper] * rg1[:courierzipper] * H * rg1[:courierzipper]' * rg2[:courierzipper]' * rg3[:courierzipper]' * rg4[:emptyzipper]' |> crystalspectrum |> visualize

function momentumoccupations(correlations::FockMap)::FockMap
    kcorrelations::Base.Generator = correlations |> crystalsubmaps
    crystal::Crystal = correlations |> getoutspace |> getcrystal
    center::Offset = crystal |> getspace |> getorigin
    tracecrystal::Crystal = Crystal(center |> Subset, crystal |> size)
    mode::Mode = Mode(:pos => center)

    function tracing(k::Momentum, corr::FockMap)::FockMap
        space::FockSpace = mode |> setattr(:offset => k) |> FockSpace
        return FockMap(space, space, [corr |> tr][:, :] |> SparseMatrixCSC) / dimension(corr |> getoutspace)
    end

    occupations::FockMap = directsum(tracing(k, corr) for (k, corr) in kcorrelations)
    fockspace::FockSpace = FockSpace(occupations |> getoutspace, reflected=tracecrystal)
    return FockMap(occupations, inspace=fockspace, outspace=fockspace)
end

using Compat
function Zipper.visualize(spectrum::CrystalSpectrum{2}; title="", toppadding::Bool = true)
    kspectrum::Dict{Momentum} = Dict(k => ([spectrum.eigenvalues[m] for m in modes] |> sort) for (k, modes) in spectrum.eigenmodes)
    mesh::Matrix{Momentum} = spectrum.crystal |> brillouinmesh
    plottingdata::Matrix{Vector} = map(k -> haskey(kspectrum, k) ? kspectrum[k] : [], mesh)
    bandcount::Integer = map(v -> v |> length, plottingdata) |> maximum
    function padding(v::Vector)::Vector
        return toppadding ? vcat(v, repeat([NaN], bandcount - length(v))) : vcat(repeat([NaN], bandcount - length(v)), v)
    end
    paddeddata::Matrix{Vector} = map(v -> v |> padding, plottingdata)
    plottingspectrum::Array = paddeddata |> stack
    layout::Layout = Layout(title=title)
    plot([contour(z=plottingspectrum[n, :, :]) for n in axes(plottingspectrum, 1)], layout)
end

filled = rg1[:filledzipper]' * rg1[:filledzipper]
empty = rg1[:emptyzipper]' * rg1[:emptyzipper]

filled1 = rg1[:courierzipper]' * rg2[:filledzipper]' * rg2[:filledzipper] * rg1[:courierzipper]
empty1 = rg1[:courierzipper]' * rg2[:emptyzipper]' * rg2[:emptyzipper] * rg1[:courierzipper]

filled2 = rg1[:courierzipper]' * rg2[:courierzipper]' * rg3[:filledzipper]' * rg3[:filledzipper] * rg2[:courierzipper] * rg1[:courierzipper]
empty2 = rg1[:courierzipper]' * rg2[:courierzipper]' * rg3[:emptyzipper]' * rg3[:emptyzipper] * rg2[:courierzipper] * rg1[:courierzipper]

filled3 = rg1[:courierzipper]' * rg2[:courierzipper]' * rg3[:courierzipper]' * rg4[:filledzipper]' * rg4[:filledzipper] * rg3[:courierzipper] * rg2[:courierzipper] * rg1[:courierzipper]
empty3 = rg1[:courierzipper]' * rg2[:courierzipper]' * rg3[:courierzipper]' * rg4[:emptyzipper]' * rg4[:emptyzipper] * rg3[:courierzipper] * rg2[:courierzipper] * rg1[:courierzipper]

core = rg1[:courierzipper]' * rg2[:courierzipper]' * rg3[:courierzipper]' * rg4[:courierzipper]' * rg4[:courierzipper] * rg3[:courierzipper] * rg2[:courierzipper] * rg1[:courierzipper]

occ1 = momentumoccupations(filled) + momentumoccupations(empty)
occ2 = momentumoccupations(filled1) + momentumoccupations(empty1)
occ3 = momentumoccupations(filled2) + momentumoccupations(empty2)
occ4 = momentumoccupations(filled3) + momentumoccupations(empty3)
occcore = momentumoccupations(core)

occ1 |> crystalspectrum |> visualize
occ2 |> crystalspectrum |> visualize
occ3 |> crystalspectrum |> visualize
occ4 |> crystalspectrum |> visualize
occcore |> crystalspectrum |> visualize

occ1 + occ2 + occ3 + occ4 + occcore |> crystalspectrum |> visualize

occ1 + occ2 + occ3 + occcore |> crystalspectrum |> visualize

frozenocc = momentumoccupations(rg1[:filledcorrelations] + rg1[:emptycorrelations]) |> crystalspectrum
frozenocc |> FockMap |> eigspech |> visualize

occ0 = momentumoccupations(C)
occ1 = momentumoccupations(rg1[:courierzipper]' * rg1[:correlations] * rg1[:courierzipper])
occ2 = momentumoccupations(rg1[:courierzipper]' * rg2[:courierzipper]' * rg2[:correlations] * rg2[:courierzipper] * rg1[:courierzipper]) - occ1
occ3 = momentumoccupations(rg1[:courierzipper]' * rg2[:courierzipper]' * rg3[:courierzipper]' * rg3[:correlations] * rg3[:courierzipper] * rg2[:courierzipper] * rg1[:courierzipper]) - occ2 - occ1

occ0 |> crystalspectrum |> visualize
occ1 |> crystalspectrum |> visualize
occ2 |> crystalspectrum |> visualize
occ3 |> crystalspectrum |> visualize

occ1 + occ2 + occ3 |> crystalspectrum |> visualize

"""

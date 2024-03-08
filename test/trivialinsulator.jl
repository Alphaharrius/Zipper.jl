using Plotly, SmithNormalForm, LinearAlgebra, OrderedCollections, SparseArrays, Combinatorics
using Zipper

fiodir("/Users/alphaharrius/ZERData")
setmaxthreads(Threads.nthreads())

triangular = RealSpace([sqrt(3)/2 -1/2; 0. 1.]')
kspace = convert(MomentumSpace, triangular)

pa = [1/3, 2/3] ∈ triangular
pb = [2/3, 1/3] ∈ triangular
pc = (pa + pb) / 2
spatialsnappingcalibration((pa, pb, pc))

c6 = pointgrouptransform([cos(π/3) -sin(π/3); sin(π/3) cos(π/3)])

unitcell = Subset(pa, pb)
crystal = Crystal(unitcell, [24, 24])
reciprocalhashcalibration(crystal.sizes)

modes::Subset{Mode} = quantize(unitcell, 1)|>orderedmodes
m0, m1 = members(modes)

tₙ = ComplexF64(-1)
t_a = ComplexF64(-0.1)
t_b = ComplexF64(-0.9)

bonds::FockMap = bondmap([
    (m0, m0) => t_a,
    (m1, m1) => t_b,
    (m0, m1) => tₙ,
    (m0, setattr(m1, :r => Point([-1, 0], triangular))) => tₙ,
    (m0, setattr(m1, :r => Point([0, 1], triangular))) => tₙ])

energyspectrum = @time computeenergyspectrum(bonds, crystal=crystal)
energyspectrum|>visualize

groundstates::CrystalSpectrum = groundstatespectrum(energyspectrum, perunitcellfillings=1)
groundstates|>visualize

groundstateprojector = groundstates|>crystalprojector
correlations = idmap(groundstateprojector|>getoutspace) - groundstateprojector

H = CrystalFockMap(energyspectrum)

function zer(correlations)
    @info("Starting RG...")
    crystalfock = correlations|>getoutspace

    scale = Scale([2 0; 0 2], crystalfock|>getcrystal|>getspace)
    @info("Performing blocking...")
    @info("Generating blocking transformation...")
    blocker = @time scale * crystalfock
    @info("Performing blocking on correlations...")
    blockedcorrelations = @time blocker * correlations * blocker'
    blockedcrystalfock = blockedcorrelations|>getoutspace
    blockedcrystal::Crystal = blockedcrystalfock|>getcrystal
    blockedspace::RealSpace = blockedcrystal|>getspace

    distillregion::Region = getsphericalregion(crystal=blockedcrystal, radius=1, metricspace=blockedspace|>orthospace)
    frozenseedingfock::RegionFock = quantize(distillregion, 1)

    @info("Computing global distill Hamiltonian...")
    globaldistiller = @time globaldistillerhamiltonian(
        correlations=blockedcorrelations,
        restrictspace=frozenseedingfock,
        localisometryselectionstrategy=frozenselectionbycount(3))

    globaldistillerspectrum = globaldistiller|>crystalspectrum

    @info("Performing groupbands...")
    distillresult = @time groupbands(globaldistillerspectrum, :courier => v -> abs(v) < 1e-5, :empty => v -> v > 1e-5, :filled => v -> v < -1e-5)

    couriersamplingregion::Region = getsphericalregion(crystal=blockedcrystal, radius=2, metricspace=blockedcrystal|>getspace|>orthospace)
    courierseedingcenter::Offset = [2/3, 1/3] ∈ blockedspace
    courierseedingregion::Region = couriersamplingregion|>filter(p -> norm((blockedspace|>orthospace) * (p - courierseedingcenter)) < 0.9)
    courierseedingfock::RegionFock = quantize(courierseedingregion, 1)

    c3 = c6^2 |> recenter(courierseedingcenter)

    function renormalizecourier()
        blockedcourierprojector = distillresult[:courier]|>crystalprojector
        blockedcouriercorrelation = idmap(blockedcourierprojector|>getoutspace) - blockedcourierprojector

        localcourierseed = findlocalspstates(
            statecorrelations=blockedcouriercorrelation,
            regionfock=courierseedingfock,
            symmetry=c3,
            spectrumextractpredicate=v -> v < 5e-2,
            statecrystalfock=blockedcrystalfock,
            linearindependencethreshold=0.0002)[1]
        fullcourierseed = localcourierseed + (c6 * localcourierseed.outspace) * localcourierseed * (c6 * localcourierseed.inspace)'

        crystalcourierseeds = crystalisometries(localisometry=fullcourierseed, crystalfock=blockedcrystalfock, addinspacemomentuminfo=true)
        wanniercourierisometry = wannierprojection(
            crystalisometries=distillresult[:courier].eigenvectors, crystal=blockedcrystal, crystalseeds=crystalcourierseeds|>Dict)

        couriercorrelations = wanniercourierisometry' * blockedcorrelations * wanniercourierisometry
        couriercorrelationspectrum = couriercorrelations |> crystalspectrum
        nonpurifiedcorrelationspectrum = couriercorrelationspectrum

        purifiedcorrelationspectrum = couriercorrelationspectrum |> roundingpurification
        couriercorrelations = purifiedcorrelationspectrum |> CrystalFockMap

        return couriercorrelations, wanniercourierisometry, nonpurifiedcorrelationspectrum
    end
    
    function renormalizefilled()
        blockedfilledprojector = distillresult[:filled]|>crystalprojector
        blockedfilledcorrelations = idmap(blockedfilledprojector|>getoutspace) - blockedfilledprojector

        filledseed = findlocalspstates(
            statecorrelations=blockedfilledcorrelations,
            regionfock=frozenseedingfock,
            symmetry=c6,
            spectrumextractpredicate=v -> v < 1e-2,
            degeneracythreshold=1e-3,
            statecrystalfock=blockedcrystalfock,
            linearindependencethreshold=0.0002)[3]

        crystalfilledseeds = crystalisometries(localisometry=filledseed, crystalfock=blockedcrystalfock, addinspacemomentuminfo=true)
        wannierfilledisometry = wannierprojection(
            crystalisometries=distillresult[:filled].eigenvectors, crystal=blockedcrystal, crystalseeds=crystalfilledseeds|>Dict)

        filledcorrelations = wannierfilledisometry' * blockedcorrelations * wannierfilledisometry
        return filledcorrelations, wannierfilledisometry
    end

    function renormalizeempty()
        blockedemptyprojector = distillresult[:empty]|>crystalprojector
        blockedemptycorrelations = idmap(blockedemptyprojector|>getoutspace) - blockedemptyprojector

        emptyseed = findlocalspstates(
            statecorrelations=blockedemptycorrelations,
            regionfock=frozenseedingfock,
            symmetry=c6,
            spectrumextractpredicate=v -> v < 1e-2,
            degeneracythreshold=1e-3,
            statecrystalfock=blockedcrystalfock,
            linearindependencethreshold=0.0002)[3]

        crystalemptyseeds = crystalisometries(localisometry=emptyseed, crystalfock=blockedcrystalfock, addinspacemomentuminfo=true)

        wannieremptyisometry = wannierprojection(
            crystalisometries=distillresult[:empty].eigenvectors, crystal=blockedcrystal, crystalseeds=crystalemptyseeds|>Dict)

        emptycorrelations = wannieremptyisometry' * blockedcorrelations * wannieremptyisometry

        return emptycorrelations, wannieremptyisometry
    end

    @info ("Renormalizing courier...")
    couriercorrelations, wanniercourierisometry, nonpurifiedcorrelationspectrum = @time renormalizecourier()
    @info ("Renormalizing filled...")
    filledcorrelations, wannierfilledisometry = @time renormalizefilled()
    @info ("Renormalizing empty...")
    emptycorrelations, wannieremptyisometry = @time renormalizeempty()

    return Dict(
        :blocker => blocker,
        :correlations => couriercorrelations,
        :emptyisometry => wannieremptyisometry,
        :filledisometry => wannierfilledisometry,
        :courierisometry => wanniercourierisometry,
        :globaldistiller => globaldistiller,
        :filledcorrelations => filledcorrelations, 
        :emptycorrelations => emptycorrelations,
        :frozenseedingfock => frozenseedingfock,
        :nonpurifiedcorrelationspectrum => nonpurifiedcorrelationspectrum)
end

rg1 = @time zer(correlations)

rg2 = @time zer(rg1[:correlations])

# rg3 = @time zer(rg2[:correlations])

# rg4 = @time zer(rg3[:correlations])

# rg5 = @time zer(rg4[:correlations])

# rg6 = @time zer(rg5[:correlations])

# rg7 = zer(rg6[:correlations])

c6g = c6 * getoutspace(rg4[:globaldistiller])
rg4[:globaldistiller] + c6g'rg4[:globaldistiller]*c6g|>crystalspectrum|>visualize

H = CrystalFockMap(energyspectrum)

function getphysspreaddata(H)
    restrictregion = getsphericalregion(crystal=crystal, radius=6, metricspace=euclidean(RealSpace, 2))
    refpos = [2/3, 1/3] ∈ getspace(crystal)
    centerregion = Subset(r for r in restrictregion if ((refpos-r)|>euclidean|>norm) <= 3.5)
    restrictfock = quantize(centerregion, 1)
    Ft = fourier(H|>getoutspace, restrictfock) / (H|>getoutspace|>getcrystal|>vol|>sqrt)
    rgHrs = Ft' * H * Ft
    visualmode = quantize(refpos|>Subset, 1)|>first
    visualH = rgHrs[:, visualmode]
    data = [((getpos(visualmode)-getpos(m)|>euclidean|>norm), (visualH[m, :]|>rep)[1, 1]|>abs) for m in visualH|>getoutspace]
    return sort(data, by=first)
end

function getspreaddata(rgdata, H)
    courierzipper = rgdata[:blocker]' * rgdata[:courierisometry]
    rgH = courierzipper' * H * courierzipper
    rgcrystal = rgH|>getoutspace|>getcrystal
    ospace = H|>getoutspace|>getcrystal|>getspace|>orthospace
    restrictregion = getsphericalregion(crystal=rgcrystal, radius=12, metricspace=ospace)
    refpos = [2/3, 1/3] ∈ getspace(rgcrystal)
    centerregion = Subset(r for r in restrictregion if (ospace*(refpos-r)|>norm) < 7)
    restrictfock = quantize(centerregion, 1)
    Ft = fourier(rgH|>getoutspace, restrictfock) / (rgH|>getoutspace|>getcrystal|>vol|>sqrt)
    rgHrs = Ft' * rgH * Ft
    visualmode = quantize(refpos|>Subset, 1)|>first
    rgvisualH = rgHrs[:, visualmode]
    data = [((getpos(visualmode)-getpos(m)|>euclidean|>norm), (rgvisualH[m, :]|>rep)[1, 1]|>abs) for m in rgvisualH|>getoutspace]
    return sort(data, by=first), rgH
end

data0 = getphysspreaddata(H)
data1, rg1H = getspreaddata(rg1, H)
data2, rg2H = getspreaddata(rg2, rg1H)
data3, rg3H = getspreaddata(rg3, rg2H)
data4, rg4H = getspreaddata(rg4, rg3H)
data5, rg5H = getspreaddata(rg5, rg4H)
data6, rg6H = getspreaddata(rg6, rg5H)

pdata0 = data0[2:end]
pdata1 = data1[2:end]
pdata2 = data2[2:end]
pdata3 = data3[2:end]
pdata4 = data4[2:end]
pdata5 = data5[2:end]
pdata6 = data6[2:end]

rg2H = fioload("rg2H")
rg2H|>crystalspectrum|>visualize

pdata0 = fioload("pdata0")
pdata1 = fioload("pdata1")
pdata2 = fioload("pdata2")
pdata3 = fioload("pdata3")
pdata4 = fioload("pdata4")
pdata5 = fioload("pdata5")

marker0 = attr(symbol="circle-open", size=12, line_width=2, color="Red")
marker1 = attr(symbol="circle-open", size=12, line_width=2, color="OrangeRed")
marker2 = attr(symbol="circle-open", size=12, line_width=2, color="Gold")
marker3 = attr(symbol="circle-open", size=12, line_width=2, color="LimeGreen")
marker4 = attr(symbol="circle-open", size=12, line_width=2, color="DeepSkyBlue")
marker5 = attr(symbol="circle-open", size=12, line_width=2, color="Blue")

plot([
    scatter(x=[d[1] for d in pdata0], y=[d[2] for d in pdata0], mode="markers", opacity=0.1, marker=marker, name="Physical"),
    scatter(x=[d[1] for d in pdata0], y=[d[2] for d in pdata1], mode="markers", opacity=0.2, marker=marker, name="RG1"),
    scatter(x=[d[1] for d in pdata0], y=[d[2] for d in pdata2], mode="markers", opacity=0.3, marker=marker, name="RG2"),
    scatter(x=[d[1] for d in pdata0], y=[d[2] for d in pdata3], mode="markers", opacity=0.4, marker=marker, name="RG3"),
    scatter(x=[d[1] for d in pdata0], y=[d[2] for d in pdata4], mode="markers", opacity=0.6, marker=marker, name="RG4"),
    scatter(x=[d[1] for d in pdata0], y=[d[2] for d in pdata5], mode="markers", opacity=0.8, marker=marker, name="RG5")])

plot([
    scatter(x=[d[1]|>log for d in pdata0], y=[d[2]|>log for d in pdata0], mode="markers", marker=marker5, name="Physical"),
    scatter(x=[d[1]|>log for d in pdata0], y=[d[2]|>log for d in pdata1], mode="markers", marker=marker4, name="RG1"),
    scatter(x=[d[1]|>log for d in pdata0], y=[d[2]|>log for d in pdata2], mode="markers", marker=marker3, name="RG2"),
    scatter(x=[d[1]|>log for d in pdata0], y=[d[2]|>log for d in pdata3], mode="markers", marker=marker2, name="RG3"),
    scatter(x=[d[1]|>log for d in pdata0], y=[d[2]|>log for d in pdata4], mode="markers", marker=marker1, name="RG4"),
    scatter(x=[d[1]|>log for d in pdata0], y=[d[2]|>log for d in pdata5], mode="markers", marker=marker0, name="RG5")],
    Layout(yaxis_range=[-10, 1]))

plot(scatter(y=[d[2] for d in data1]))

visualize(rg3Hrs[:, 135]|>RegionState, markersizemultiplier=20, markersizescaling=0.1)

lineregionfock = RegionFock(m for m in restrictfock if Integer((m|>getattr(:r)|>vec)[2]|>round) == 1)
lineregionfock|>getregion|>visualize

rg1Hline = rg1Hrs[lineregionfock, lineregionfock]
anchormode = sort([m=>m|>getpos|>norm for m in lineregionfock], by=last)|>last|>first
datasource = rg1Hline[:, anchormode]
[(getpos(anchormode)-getpos(m)|>euclidean|>norm, (datasource[m, :]|>rep)[1, 1]|>real) for m in datasource|>getoutspace]
rg1Hline|>visualize
rg1Hline[anchormode, anchormode]|>rep

entanglemententropy(rg1[:filledcorrelations]|>crystalspectrum) / (rg1[:filledcorrelations]|>getoutspace|>getcrystal|>vol)
entanglemententropy(rg2[:filledcorrelations]|>crystalspectrum) / (rg2[:filledcorrelations]|>getoutspace|>getcrystal|>vol)
entanglemententropy(rg3[:filledcorrelations]|>crystalspectrum) / (rg3[:filledcorrelations]|>getoutspace|>getcrystal|>vol)
entanglemententropy(rg4[:filledcorrelations]|>crystalspectrum) / (rg4[:filledcorrelations]|>getoutspace|>getcrystal|>vol)
entanglemententropy(groundstatespectrum(rg4[:nonpurifiedcorrelationspectrum], perunitcellfillings=1)) / (rg4[:correlations]|>getoutspace|>getcrystal|>vol)

rg5[:globaldistiller]|>crystalspectrum|>visualize

rg1[:entanglemententropy] / (rg1[:correlations]|>getoutspace|>getcrystal|>vol)
rg2[:entanglemententropy] / (rg2[:correlations]|>getoutspace|>getcrystal|>vol)
rg3[:entanglemententropy] / (rg3[:correlations]|>getoutspace|>getcrystal|>vol)
rg4[:entanglemententropy] / (rg4[:correlations]|>getoutspace|>getcrystal|>vol)

rg1filledC = rg1[:filledzipper]' * rg1[:filledzipper]
rg2filledC = rg1[:courierzipper]' * rg2[:filledzipper]' * rg2[:filledzipper] * rg1[:courierzipper]
rg3filledC = rg1[:courierzipper]' * rg2[:courierzipper]' * rg3[:filledzipper]' * rg3[:filledzipper] * rg2[:courierzipper] * rg1[:courierzipper]
rg4filledC = rg1[:courierzipper]' * rg2[:courierzipper]' * rg3[:courierzipper]' * rg4[:filledzipper]' * rg4[:filledzipper] * rg3[:courierzipper] * rg2[:courierzipper] * rg1[:courierzipper]
rg4P = groundstatespectrum(rg4[:correlations]|>crystalspectrum, perunitcellfillings=1)|>crystalprojector
rg4C = rg1[:courierzipper]' * rg2[:courierzipper]' * rg3[:courierzipper]' * rg4[:courierzipper]' * rg4P * rg4[:courierzipper] * rg3[:courierzipper] * rg2[:courierzipper] * rg1[:courierzipper]

rgP = rg1filledC + rg2filledC + rg3filledC + rg4filledC + rg4C
rgC = idmap(rgP|>getoutspace) - rgP
visualize(rgC, colrange=1:512, rowrange=1:512)
rgC|>crystalspectrum|>visualize

visualize(rgC - C, colrange=1:512, rowrange=1:512)
M = rgC - C
(M * M' |> tr) / 

function circularregionmodes(origin::Offset, physicalmodes::Subset{Mode}, radius::Number)::Subset{Mode}
    currentspace::RealSpace = blockedcrystal |> getspace |> orthospace
    physicalnorm = m -> lineartransform(currentspace, m |> getpos) |> norm
    return filter(p -> physicalnorm(p - origin) < radius, physicalmodes)
end

# blockedcrystal = rg1[:correlations]|>getoutspace|>getcrystal
# crystalpoints::Subset{Offset} = latticepoints(blockedcrystal)
# samplepoints::Subset{Offset} = crystalpoints + c6^2 * crystalpoints + c6^4 * crystalpoints
# blockedmodes::Subset{Mode} = quantize(:b, blockedcrystal.unitcell, 1)
# physicalmodes::Subset{Mode} = spanoffset(blockedmodes, samplepoints)


# frozenseedingmodes::Subset{Mode} = circularregionmodes(triangular |> getorigin, physicalmodes, 8)
# frozenseedingfock::FockSpace = RegionFock(frozenseedingmodes)


# regionalrestriction(rg1[:filledisometry], frozenseedingfock) |> visualize

# rg1[:courierzipper]' * rg2[:blocker]' * rg2[:filledisometry]

# regionalrestriction(rg1[:courierisometry] * rg2[:blocker]' * rg2[:filledisometry], frozenseedingfock) |> visualize
# rg1[:courierisometry] * rg2[:courierzipper]' * rg3[:blocker]' * rg3[:filledisometry]
# regionalrestriction(rg1[:courierisometry] * rg2[:courierzipper]' * rg3[:blocker]' * rg3[:filledisometry], frozenseedingfock) |> visualize

# rg3[:correlations] |> crystalspectrum |> visualize

# regionalrestriction(rg3[:blocker]' * rg3[:filledisometry], rg2[:frozenseedingfock]) |> visualize

# rg1[:entanglemententropy]

# rg1[:globaldistiller] |> crystalspectrum |> visualize
# rg1[:filledisometry] |> getinspace |> unitcellfock |> modeattrs

# rg1[:filledcorrelations] |> crystalspectrum |> visualize

H = energyspectrum |> FockMap

function Chern_number_multiband(state::CrystalSpectrum)
    statevectors::Dict{Momentum, FockMap} = state |> geteigenvectors
    statecrystal = state |> getcrystal
    Δkx = (1 / size(statecrystal)[1], 0) ∈ kspace
    Δky = (0, 1 / size(statecrystal)[2]) ∈ kspace
    Δk = [Δkx, Δky]
    
    function LinkVariable(k::Momentum, direction)
        delta_k = Δk[direction]
        overlap = (statevectors[k] |> rep)' * (statevectors[k+delta_k|> basispoint]|> rep)
        U_mu = det(overlap)
        U_mu /= abs(U_mu)
        return U_mu
    end
    
    function Berry_phase(k::Momentum)
        return LinkVariable(k, 1) * LinkVariable(k + Δk[1]|> basispoint, 2) / (LinkVariable(k + Δk[2]|> basispoint, 1) * LinkVariable(k, 2))
    end

    function Berry_curvature(k::Momentum)
        return -1im*log(Berry_phase(k))
    end

    Chern = 0 
    for kpoints in statevectors
        k = kpoints[1]
        Chern += Berry_curvature(k)
    end
    return Chern/(2*pi)
end

function berrycurvaturemultiband(state::CrystalSpectrum)
    statevectors::Dict{Momentum, FockMap} = state |> geteigenvectors
    statecrystal = state |> getcrystal
    Δkx = (1 / size(statecrystal)[1], 0) ∈ kspace
    Δky = (0, 1 / size(statecrystal)[2]) ∈ kspace
    Δk = [Δkx, Δky]
    
    function LinkVariable(k::Momentum, direction)
        delta_k = Δk[direction]
        overlap = (statevectors[k] |> rep)' * (statevectors[k+delta_k|> basispoint]|> rep)
        U_mu = det(overlap)
        U_mu /= abs(U_mu)
        return U_mu
    end
    
    function Berry_phase(k::Momentum)
        return LinkVariable(k, 1) * LinkVariable(k + Δk[1]|> basispoint, 2) / (LinkVariable(k + Δk[2]|> basispoint, 1) * LinkVariable(k, 2))
    end

    function Berry_curvature(k::Momentum)
        return -1im*log(Berry_phase(k))
    end

    function kberrycurvature(k)
        mode = Mode(:k => k, :b => statecrystal|>getspace|>getorigin)
        fockspace = FockSpace(mode)
        return FockMap(fockspace, fockspace, [Berry_curvature(k)][:, :])
    end

    fockmap = directsum(kberrycurvature(k) for (k, _) in statevectors)
    barecrystal = Crystal(statecrystal|>getspace|>getorigin, statecrystal|>size)
    crystalfock = FockSpace(fockmap|>getoutspace, reflected=barecrystal)
    return FockMap(fockmap, inspace=crystalfock, outspace=crystalfock)
end
:b
scale

Chern_number_multiband(groundstates)
berrycurvaturemap = berrycurvaturemultiband(groundstates)
visualize(berrycurvaturemap|>crystalspectrum, usecontour=true, title="berrycurvature-conventional-GS")

rg1courierprojector = rg1[:courierzipper]' * rg1[:courierzipper]
rg1H = rg1courierprojector * H * rg1courierprojector'
rg1gs = groundstatespectrum(rg1H |> crystalspectrum, perunitcellfillings=1)
visualize(rg1gs, usecontour=true, title="gsH-RG1")

rg1H = rg1[:courierzipper] * H * rg1[:courierzipper]'
rg1gs = groundstatespectrum(rg1H|>crystalspectrum, perunitcellfillings=1)
visualize(rg1gs, usecontour=true, title="gsH-RG1-RG1space")

using SparseArrays
function quantummetricdeterminants(state::CrystalSpectrum)
    spatialorigin = state|>getcrystal|>getspace|>getorigin
    barecrystal = Crystal(spatialorigin|>Subset, state|>getcrystal|>size)
    
    function makedetblock(k, kgeometrictensor)
        mode = Mode(:k => k, :b => spatialorigin)
        fockspace = FockSpace(mode)
        return FockMap(fockspace, fockspace, [kgeometrictensor|>rep|>imag|>det][:, :]|>SparseMatrixCSC)
    end

    kqgtensors = computequantumgeometrictensor2d(state)
    product = directsum(makedetblock(k, kgeometrictensor) for (k, kgeometrictensor) in kqgtensors)
    fockspace = FockSpace(product|>getoutspace, reflected=barecrystal)
    return FockMap(product, inspace=fockspace, outspace=fockspace)
end
:b
rg1det = quantummetricdeterminants(rg1gs)
visualize(rg1det, colrange=512:1024, rowrange=512:1024)
visualize(rg1det|>crystalspectrum, usecontour=true, title="qgtdet-RG1-RG1space")

rg1occ = momentumoccupations(rg1[:correlations])
rg1occ|>crystalspectrum|>geteigenvalues

rg1qm = computequantummetric(rg1gs)
visualize(rg1qm[:gxx], colrange=512:1024, rowrange=512:1024)

visualize(rg1qm[:gxx]|>crystalspectrum, usecontour=true, title="gxx-RG1-RG1space")
visualize(rg1qm[:gyy]|>crystalspectrum, usecontour=true, title="gyy-RG1-RG1space")
visualize(rg1qm[:gxy]|>crystalspectrum, usecontour=true, title="gxy-RG1-RG1space")
visualize(rg1qm[:gyx]|>crystalspectrum, usecontour=true, title="gyx-RG1-RG1space")

rg2courierisometry = rg2[:courierzipper] * rg1[:courierzipper]
rg2courierprojector = rg2courierisometry' * rg2courierisometry
rg2H = rg2courierprojector * H * rg2courierprojector'
rg2gs = groundstatespectrum(rg2H |> crystalspectrum, perunitcellfillings=1)
visualize(rg2gs, usecontour=true, title="gsH-RG2")

rg2H = rg2[:courierzipper] * rg1H * rg2[:courierzipper]'
rg2gs = groundstatespectrum(rg2H|>crystalspectrum, perunitcellfillings=1)
visualize(rg2gs, usecontour=true, title="gsH-RG2-RG2space")

rg2det = quantummetricdeterminants(rg2gs)
visualize(rg2det|>crystalspectrum, usecontour=true, title="qgtdet-RG2-RG2space")

rg2qm = computequantummetric(rg2gs)
visualize(rg2qm[:gxx]|>crystalspectrum, usecontour=true, title="gxx-RG2-RG2space")
visualize(rg2qm[:gyy]|>crystalspectrum, usecontour=true, title="gyy-RG2-RG2space")

rg3courierisometry = rg3[:courierzipper] * rg2[:courierzipper] * rg1[:courierzipper]
rg3courierprojector = rg3courierisometry' * rg3courierisometry
rg3H = rg3courierprojector * H * rg3courierprojector'
rg3gs = groundstatespectrum(rg3H |> crystalspectrum, perunitcellfillings=1)
visualize(rg3gs, usecontour=true, title="gsH-RG3")

rg3H = rg3[:courierzipper] * rg2H * rg3[:courierzipper]'
rg3gs = groundstatespectrum(rg3H|>crystalspectrum, perunitcellfillings=1)
visualize(rg3gs, usecontour=true, title="gsH-RG3-RG3space")

rg3det = quantummetricdeterminants(rg3gs)
visualize(rg3det|>crystalspectrum, usecontour=true, title="qgtdet-RG3-RG3space")

rg3qm = computequantummetric(rg3gs)
visualize(rg3qm[:gxx]|>crystalspectrum, usecontour=true, title="gxx-RG3-RG3space")
visualize(rg3qm[:gyy]|>crystalspectrum, usecontour=true, title="gyy-RG3-RG3space")

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

visualize(occ1|>crystalspectrum, usecontour=true, title="occ1")
visualize(occ2|>crystalspectrum, usecontour=true, title="occ2")
visualize(occ3|>crystalspectrum, usecontour=true, title="occ3")
visualize(occ4|>crystalspectrum, usecontour=true, title="occ4")
visualize(occcore|>crystalspectrum, usecontour=true, title="occcore")

occ1 + occ2 + occ3 + occ4 + occcore |> crystalspectrum |> visualize

occ1 + occ2 + occ3 + occcore |> crystalspectrum |> visualize

frozenocc = momentumoccupations(rg1[:filledcorrelations] + rg1[:emptycorrelations]) |> crystalspectrum
frozenocc |> FockMap |> eigspech |> visualize

occ0 = momentumoccupations(C)
occ1 = momentumoccupations(rg1[:courierzipper]' * rg1[:correlations] * rg1[:courierzipper])
occ2 = momentumoccupations(rg1[:courierzipper]' * rg2[:courierzipper]' * rg2[:correlations] * rg2[:courierzipper] * rg1[:courierzipper]) - occ1
occ3 = momentumoccupations(rg1[:courierzipper]' * rg2[:courierzipper]' * rg3[:courierzipper]' * rg3[:correlations] * rg3[:courierzipper] * rg2[:courierzipper] * rg1[:courierzipper]) - occ2 - occ1

rg1C = rg1[:courierzipper]' * rg1[:correlations] * rg1[:courierzipper]
rg2C = rg1[:courierzipper]' * rg2[:courierzipper]' * rg2[:correlations] * rg2[:courierzipper] * rg1[:courierzipper]
rg3C = rg1[:courierzipper]' * rg2[:courierzipper]' * rg3[:courierzipper]' * rg3[:correlations] * rg3[:courierzipper] * rg2[:courierzipper] * rg1[:courierzipper]

visualize(rg1C, colrange=1:512, rowrange=1:512)
visualize(rg2C, colrange=1:512, rowrange=1:512)
visualize(rg3C, colrange=1:512, rowrange=1:512)

occ0 |> crystalspectrum |> visualize
occ1 |> crystalspectrum |> visualize
occ2 |> crystalspectrum |> visualize
occ3 |> crystalspectrum |> visualize

occ1 + occ2 + occ3 |> crystalspectrum |> visualize

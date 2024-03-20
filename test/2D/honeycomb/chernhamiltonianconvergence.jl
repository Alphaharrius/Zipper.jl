using LinearAlgebra
using Zipper, Plots
plotlyjs()

setmaxthreads(Threads.nthreads())

@info "Preparing environment..."
triangular = RealSpace([sqrt(3)/2 -1/2; 0. 1.]')
kspace = convert(MomentumSpace, triangular)

pa = [1/3, 2/3] ∈ triangular
pb = [2/3, 1/3] ∈ triangular
pc = (pa + pb) / 2
spatialsnappingcalibration((pa, pb, pc))

c6 = pointgrouptransform([cos(π/3) -sin(π/3); sin(π/3) cos(π/3)])
c3 = c6^2

unitcell = Subset(pa, pb)
crystal = Crystal(unitcell, [192, 192])
reciprocalhashcalibration(crystal.sizes)

@info "Loading data..."
fiodir("/Users/alphaharrius/ZERData/chern192x192")
hamiltonian = fioload("hamiltonian")
correlations = fioload("correlations")

fiodir("/Users/alphaharrius/ZERData/chern192x192/rg1")
rg1correlations = fioload("couriercorrelations")
rg1block = fioload("block")
rg1zipper = fioload("wanniercourierisometry")' * rg1block

fiodir("/Users/alphaharrius/ZERData/chern192x192/rg2")
rg2correlations = fioload("couriercorrelations")
rg2block = fioload("block")
rg2zipper = fioload("wanniercourierisometry")' * rg2block

fiodir("/Users/alphaharrius/ZERData/chern192x192/rg3")
rg3correlations = fioload("couriercorrelations")
rg3block = fioload("block")
rg3zipper = fioload("wanniercourierisometry")' * rg3block

fiodir("/Users/alphaharrius/ZERData/chern192x192/rg4")
rg4correlations = fioload("couriercorrelations")
rg4block = fioload("block")
rg4zipper = fioload("wanniercourierisometry")' * rg4block

fiodir("/Users/alphaharrius/ZERData/chern192x192/rg5")
rg5correlations = fioload("couriercorrelations")
rg5block = fioload("block")
rg5zipper = fioload("wanniercourierisometry")' * rg5block

fiodir("/Users/alphaharrius/ZERData/chern192x192/rg6")
rg6correlations = fioload("couriercorrelations")
rg6block = fioload("block")
rg6zipper = fioload("wanniercourierisometry")' * rg6block

@info "Computing physical data..."
region = getsphericalregion(crystal=crystal, radius=6, metricspace=euclidean(RealSpace, 2))
ref = getspace(crystal)*(2/3, 1/3)
centerregion = Subset(r for r in region if ((ref-r)|>euclidean|>norm) <= 2.5)
regionfock = getregionfock(hamiltonian|>getoutspace, centerregion)
transform = fourier(hamiltonian|>getoutspace, regionfock)
restricted = transform' * hamiltonian * transform / (crystal|>vol)
visualmode = getregionfock(hamiltonian|>getoutspace, ref|>Subset)|>first
r1 = getspace(crystal)*(1/3, 2/3)
m1 = getregionfock(hamiltonian|>getoutspace, r1|>Subset)|>first
h1 = restricted[visualmode, m1]|>rep
r2 = getspace(crystal)*(2/3, 1 + 1/3)
m2 = getregionfock(hamiltonian|>getoutspace, r2|>Subset)|>first
h2 = restricted[visualmode, m2]|>rep
r3 = getspace(crystal)*(1 + 1/3, 1 + 2/3)
m3 = getregionfock(hamiltonian|>getoutspace, r3|>Subset)|>first
h3 = restricted[visualmode, m3]|>rep
r4 = getspace(crystal)*(1 + 2/3, 2 + 1/3)
m4 = getregionfock(hamiltonian|>getoutspace, r4|>Subset)|>first
h4 = restricted[visualmode, m4]|>rep
r5 = getspace(crystal)*(2 + 1/3, 1 + 2/3)
m5 = getregionfock(hamiltonian|>getoutspace, r5|>Subset)|>first
h5 = restricted[visualmode, m5]|>rep
r6 = getspace(crystal)*(2 + 2/3, 2 + 1/3)
m6 = getregionfock(hamiltonian|>getoutspace, r6|>Subset)|>first
h6 = restricted[visualmode, m6]|>rep
visualize(centerregion, ref|>Subset, r1|>Subset, r2|>Subset, r3|>Subset, r4|>Subset, r5|>Subset, r6|>Subset)

@info "Computing renormalization group data..."
function computergspread(rghamiltonian)
    rgcrystal = rghamiltonian|>getoutspace|>getcrystal
    rgspace = rgcrystal|>getspace
    scale = Scale([2 0; 0 2], rgspace)
    ospace = inv(scale)*rgspace|>orthospace
    region = getsphericalregion(crystal=rgcrystal, radius=12, metricspace=ospace)
    ref = rgspace*(2/3, 1/3)
    centerregion = Subset(r for r in region if (ospace*(ref-r)|>norm) < 7)
    regionfock = getregionfock(rghamiltonian|>getoutspace, centerregion)
    transform = fourier(rghamiltonian|>getoutspace, regionfock) / (rgcrystal|>vol|>sqrt)
    restricted = transform' * rghamiltonian * transform
    visualmode = getregionfock(rghamiltonian|>getoutspace, ref|>Subset)|>first
    r1 = rgspace*(1/3, 2/3)
    m1 = getregionfock(rghamiltonian|>getoutspace, r1|>Subset)|>first
    h1 = restricted[visualmode, m1]|>rep
    r2 = rgspace*(2/3, 1 + 1/3)
    m2 = getregionfock(rghamiltonian|>getoutspace, r2|>Subset)|>first
    h2 = restricted[visualmode, m2]|>rep
    r3 = rgspace*(1 + 1/3, 1 + 2/3)
    m3 = getregionfock(rghamiltonian|>getoutspace, r3|>Subset)|>first
    h3 = restricted[visualmode, m3]|>rep
    r4 = rgspace*(1 + 2/3, 2 + 1/3)
    m4 = getregionfock(rghamiltonian|>getoutspace, r4|>Subset)|>first
    h4 = restricted[visualmode, m4]|>rep
    r5 = rgspace*(2 + 1/3, 1 + 2/3)
    m5 = getregionfock(rghamiltonian|>getoutspace, r5|>Subset)|>first
    h5 = restricted[visualmode, m5]|>rep
    r6 = rgspace*(2 + 2/3, 2 + 1/3)
    m6 = getregionfock(rghamiltonian|>getoutspace, r6|>Subset)|>first
    h6 = restricted[visualmode, m6]|>rep
    return Dict(:h1=>h1[1, 1], :h2=>h2[1, 1], :h3=>h3[1, 1], :h4=>h4[1, 1], :h5=>h5[1, 1], :h6=>h6[1, 1])
end

rg1hamiltonian = rg1zipper * hamiltonian * rg1zipper'
rg1spread = computergspread(rg1hamiltonian)

rg2hamiltonian = rg2zipper * rg1hamiltonian * rg2zipper'
rg2spread = computergspread(rg2hamiltonian)

rg3hamiltonian = rg3zipper * rg2hamiltonian * rg3zipper'
rg3spread = computergspread(rg3hamiltonian)

rg4hamiltonian = rg4zipper * rg3hamiltonian * rg4zipper'
rg4spread = computergspread(rg4hamiltonian)

rg5hamiltonian = rg5zipper * rg4hamiltonian * rg5zipper'
rg5spread = computergspread(rg5hamiltonian)

rg6hamiltonian = rg6zipper * rg5hamiltonian * rg6zipper'
rg6spread = computergspread(rg6hamiltonian)

band1 = [h1[1, 1], rg1spread[:h1], rg2spread[:h1], rg3spread[:h1], rg4spread[:h1], rg5spread[:h1], rg6spread[:h1]]
band2 = [h2[1, 1], rg1spread[:h2], rg2spread[:h2], rg3spread[:h2], rg4spread[:h2], rg5spread[:h2], rg6spread[:h2]]
band3 = [h3[1, 1], rg1spread[:h3], rg2spread[:h3], rg3spread[:h3], rg4spread[:h3], rg5spread[:h3], rg6spread[:h3]]
band4 = [h4[1, 1], rg1spread[:h4], rg2spread[:h4], rg3spread[:h4], rg4spread[:h4], rg5spread[:h4], rg6spread[:h4]]
band5 = [h5[1, 1], rg1spread[:h5], rg2spread[:h5], rg3spread[:h5], rg4spread[:h5], rg5spread[:h5], rg6spread[:h5]]
band6 = [h6[1, 1], rg1spread[:h6], rg2spread[:h6], rg3spread[:h6], rg4spread[:h6], rg5spread[:h6], rg6spread[:h6]]

@info "Plotting linear-scale (Re)..."
scatter(1:length(band1), band1|>real, label="$r1", legend=:bottomright)
scatter!(1:length(band2), band2|>real, label="$r2", legend=:bottomright)
scatter!(1:length(band3), band3|>real, label="$r3", legend=:bottomright)
scatter!(1:length(band4), band4|>real, label="$r4", legend=:bottomright)
scatter!(1:length(band5), band5|>real, label="$r5", legend=:bottomright)
scatter!(1:length(band6), band6|>real, label="$r6", legend=:bottomright)

@info "Plotting log-scale (Re)..."
scatter(1:length(band1), [log(v|>real|>abs) for v in band1], label="$r1", legend=:bottomright)
scatter!(1:length(band2), [log(v|>real|>abs) for v in band2], label="$r2", legend=:bottomright)
scatter!(1:length(band3), [log(v|>real|>abs) for v in band3], label="$r3", legend=:bottomright)
scatter!(1:length(band4), [log(v|>real|>abs) for v in band4], label="$r4", legend=:bottomright)
scatter!(1:length(band5), [log(v|>real|>abs) for v in band5], label="$r5", legend=:bottomright)
scatter!(1:length(band6), [log(v|>real|>abs) for v in band6], label="$r6", legend=:bottomright)

@info "Plotting linear-scale (Im)..."
scatter(1:length(band1), band1|>imag, label="$r1", legend=:bottomright)
scatter!(1:length(band2), band2|>imag, label="$r2", legend=:bottomright)
scatter!(1:length(band3), band3|>imag, label="$r3", legend=:bottomright)
scatter!(1:length(band4), band4|>imag, label="$r4", legend=:bottomright)
scatter!(1:length(band5), band5|>imag, label="$r5", legend=:bottomright)
scatter!(1:length(band6), band6|>imag, label="$r6", legend=:bottomright)
plot!(legend=false)

@info "Plotting log-scale (Im)..."
scatter(1:length(band1), [log(v|>imag|>abs) for v in band1], label="$r1", legend=:bottomright)
scatter!(1:length(band2), [log(v|>imag|>abs) for v in band2], label="$r2", legend=:bottomright)
scatter!(1:length(band3), [log(v|>imag|>abs) for v in band3], label="$r3", legend=:bottomright)
scatter!(1:length(band4), [log(v|>imag|>abs) for v in band4], label="$r4", legend=:bottomright)
scatter!(1:length(band5), [log(v|>imag|>abs) for v in band5], label="$r5", legend=:bottomright)
scatter!(1:length(band6), [log(v|>imag|>abs) for v in band6], label="$r6", legend=:bottomright)

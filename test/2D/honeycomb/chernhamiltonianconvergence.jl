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
centerregion = Subset(r for r in region if ((ref-r)|>euclidean|>norm) <= 3.5)
regionfock = getregionfock(hamiltonian|>getoutspace, centerregion)
transform = fourier(hamiltonian|>getoutspace, regionfock) / (crystal|>vol|>sqrt)
restricted = transform' * hamiltonian * transform
visualmode = getregionfock(hamiltonian|>getoutspace, ref|>Subset)|>first
visualhamiltonian = restricted[:, visualmode]
physspread = [
    ((getpos(visualmode)-getpos(m)|>euclidean|>norm), (visualhamiltonian[m, :]|>rep)[1, 1]|>abs) 
    for m in visualhamiltonian|>getoutspace]

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
    visualhamiltonian = restricted[:, visualmode]
    return [
        ((getpos(visualmode)-getpos(m)|>euclidean|>norm), (visualhamiltonian[m, :]|>rep)[1, 1]|>abs) 
        for m in visualhamiltonian|>getoutspace]
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

@info "Removing first data point..."
physspread = physspread[2:end]
rg1spread = rg1spread[2:end]
rg2spread = rg2spread[2:end]
rg3spread = rg3spread[2:end]
rg4spread = rg4spread[2:end]
rg5spread = rg5spread[2:end]
rg6spread = rg6spread[2:end]

@info "Plotting linear-scale..."
scatter(
    [x for (x, _) in physspread], [y for (_, y) in physspread], 
    label="Physical", 
    markercolor=RGBA(1, 1, 1, 0), xlabel="Distance", ylabel="Energy", legend=:topright,
    markerstrokecolor="Blue", markerstrokewidth=2, markerstrokeopacity=0.5)
scatter!(
    [x for (x, _) in physspread], [y for (_, y) in rg1spread], 
    label="RG1",
    markercolor=RGBA(1, 1, 1, 0), xlabel="Distance", ylabel="Energy", legend=:topright,
    markerstrokecolor="DeepSkyBlue", markerstrokewidth=2, markerstrokeopacity=0.5)
scatter!(
    [x for (x, _) in physspread], [y for (_, y) in rg2spread], 
    label="RG2", 
    markercolor=RGBA(1, 1, 1, 0), xlabel="Distance", ylabel="Energy", legend=:topright,
    markerstrokecolor="LimeGreen", markerstrokewidth=2, markerstrokeopacity=0.5)
scatter!(
    [x for (x, _) in physspread], [y for (_, y) in rg3spread], 
    label="RG3",
    markercolor=RGBA(1, 1, 1, 0), xlabel="Distance", ylabel="Energy", legend=:topright,
    markerstrokecolor="Gold", markerstrokewidth=2, markerstrokeopacity=0.5)
scatter!(
    [x for (x, _) in physspread], [y for (_, y) in rg4spread], 
    label="RG4", 
    markercolor=RGBA(1, 1, 1, 0), xlabel="Distance", ylabel="Energy", legend=:topright,
    markerstrokecolor="OrangeRed", markerstrokewidth=2, markerstrokeopacity=0.5)
scatter!(
    [x for (x, _) in physspread], [y for (_, y) in rg5spread], 
    label="RG5", 
    markercolor=RGBA(1, 1, 1, 0), xlabel="Distance", ylabel="Energy", legend=:topright,
    markerstrokecolor="Red", markerstrokewidth=2, markerstrokeopacity=0.5)

@info "Plotting log-scale..."
scatter(
    [x|>log for (x, _) in physspread], [y|>log for (_, y) in physspread], 
    label="Physical", 
    markercolor=RGBA(1, 1, 1, 0), xlabel="Distance", ylabel="Energy", legend=:topright,
    markerstrokecolor="Blue", markerstrokewidth=2, markerstrokeopacity=0.5)
scatter!(
    [x|>log for (x, _) in physspread], [y|>log for (_, y) in rg1spread], 
    label="RG1",
    markercolor=RGBA(1, 1, 1, 0), xlabel="Distance", ylabel="Energy", legend=:topright,
    markerstrokecolor="DeepSkyBlue", markerstrokewidth=2, markerstrokeopacity=0.5)
scatter!(
    [x|>log for (x, _) in physspread], [y|>log for (_, y) in rg2spread], 
    label="RG2", 
    markercolor=RGBA(1, 1, 1, 0), xlabel="Distance", ylabel="Energy", legend=:topright,
    markerstrokecolor="LimeGreen", markerstrokewidth=2, markerstrokeopacity=0.5)
scatter!(
    [x|>log for (x, _) in physspread], [y|>log for (_, y) in rg3spread], 
    label="RG3",
    markercolor=RGBA(1, 1, 1, 0), xlabel="Distance", ylabel="Energy", legend=:topright,
    markerstrokecolor="Gold", markerstrokewidth=2, markerstrokeopacity=0.5)
scatter!(
    [x|>log for (x, _) in physspread], [y|>log for (_, y) in rg4spread], 
    label="RG4", 
    markercolor=RGBA(1, 1, 1, 0), xlabel="Distance", ylabel="Energy", legend=:topright,
    markerstrokecolor="OrangeRed", markerstrokewidth=2, markerstrokeopacity=0.5)
scatter!(
    [x|>log for (x, _) in physspread], [y|>log for (_, y) in rg5spread], 
    label="RG5", 
    markercolor=RGBA(1, 1, 1, 0), xlabel="Distance", ylabel="Energy", legend=:topright,
    markerstrokecolor="Red", markerstrokewidth=2, markerstrokeopacity=0.5)
plot!(xlimits=(-1, 1.5), ylimits=(-10,1))

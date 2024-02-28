using Zipper

using Plots
plotlyjs()

setmaxthreads(Threads.nthreads())

fiodir("/Users/alphaharrius/ZERData/squaremetal")
inputC = fioload("inputC")
testfockmap = inputC.blocks|>first|>last
realplot = heatmap(testfockmap|>rep|>real|>Matrix)
imagplot = heatmap(testfockmap|>rep|>imag|>Matrix)
plot(realplot, imagplot, yflip=true)

function Zipper.visualize(fockmap::FockMap)
    source = fockmap|>rep
    realplot = heatmap(source|>real)
    imagplot = heatmap(source|>imag)
    return plot(realplot, imagplot, yflip=true)
end

struct PlotRegion{D}
    region::Region
end

PlotRegion(region::Region) = PlotRegion{region|>getspace|>dimension}(region)

function Zipper.visualize(region::Region; kwargs...)
    plotregion = PlotRegion(region)
    return visualize(plotregion, kwargs...)
end

function Zipper.visualize(region::PlotRegion{2}; visualspace::AffineSpace = region.region|>getspace|>euclidean)
    r = [visualspace*r|>vec for r in region.region]
    x = [v|>first for v in r]
    y = [v|>last for v in r]
    # aspect_ratio=:equal for scatter is not working as expected, it only works when setting to 1.1.
    return scatter(x, y, aspect_ratio=1.1, xlims=(minimum(x)-1, maximum(x)+1), ylims=(minimum(y)-1, maximum(y)+1))
end

crystal = inputC|>getoutspace|>getcrystal
region = getsphericalregion(crystal=crystal, radius=10, metricspace=crystal|>getspace)
visualize(region)

fiodir("/Users/alphaharrius/ZERData/squaremetal/RG1")
courierisometry = fioload("rg1courierisometry")
couriercrystal = courierisometry|>getoutspace|>getcrystal
visualregion = getsphericalregion(crystal=couriercrystal, radius=2, metricspace=couriercrystal|>getspace)
visualfock = quantize(visualregion, 1)
leftrestrict = fourier(courierisometry|>getoutspace, visualfock) / (couriercrystal|>vol|>sqrt)
rightrestrict = fourier(courierisometry|>getinspace, courierisometry|>getinspace|>unitcellfock|>RegionFock)
courierlocalstates = leftrestrict' * courierisometry * rightrestrict
courierlocalstates = RegionState(courierlocalstates)

function Zipper.visualize(
    states::RegionState{2};
    visualspace::RealSpace = states|>getoutspace|>getregion|>getspace|>euclidean,
    markersize::Integer = 1, logscale::Real = 1)

    function plotstate(order::Integer, statemap::SparseFockMap)
        r = [visualspace*getpos(mode)|>vec for mode in statemap|>getoutspace]
        x = [v|>first for v in r]
        y = [v|>last for v in r]
        f = ((statemap[mode, :]|>rep)[1, 1] for mode in statemap|>getoutspace)
        color = [convert(RGB{Float64}, HSV(angle(v) / 2π * 360, 1, 1)) for v in f]
        size = [markersize*(abs(v)^logscale) for v in f]
        hover = [trunc(v|>real, digits=3)+trunc(v|>imag, digits=3)*im|>string for v in f]
        return scatter(
            x, y, markersize=size, markercolor=color, hover=hover,
            markerstrokealpha=0.01,
            aspect_ratio=1.1, xlims=(minimum(x)-1, maximum(x)+1), ylims=(minimum(y)-1, maximum(y)+1),
            title="mode $order")
    end

    source = [(n, mode, state) for (n, (mode, state)) in states|>enumerate]
    p = plot((plotstate(n, state) for (n, _, state) in source)..., legend=false)
    @info "Showing corresponding modes..."
    showmodes(mode for (_, mode, _) in source)
    return p
end

visualize(courierlocalstates, markersize=2, logscale=0.4)

import PlotlyJS
function Zipper.visualize(spectrum::CrystalSpectrum{2}; paddirection::Symbol=:top)
    eigenmodes = spectrum|>geteigenmodes
    eigenvalues = spectrum|>geteigenvalues
    crystal = spectrum|>getcrystal
    kspectrum = Dict(
        k=>(haskey(eigenmodes, k) ? [eigenvalues[mode] for mode in eigenmodes[k]]|>sort : []) 
        for k in crystal|>brillouinzone)
    bandcount = spectrum.bandcount[2]

    function dopadding(k, values)
        pad = repeat([NaN], bandcount - length(values))
        if paddirection == :bottom
            kspectrum[k] = [pad..., values...]
        elseif paddirection == :top
            kspectrum[k] = [values..., pad...]
        else
            error("Invalid paddirection: $paddirection")
        end
    end

    if spectrum.bandcount[1] != bandcount
        watchprogress(desc="visualize(::$(spectrum|>typeof))")
        for k in spectrum|>getcrystal|>brillouinzone
            values = haskey(kspectrum, k) ? kspectrum[k] : []
            length(values) < bandcount && dopadding(k, values)
            updateprogress()
        end
        unwatchprogress()
    end
    
    getband(n) = [k=>kspectrum[k][n] for k in kspectrum|>keys]

    function plotband(band)
        K = [k|>euclidean|>vec for (k, _) in band]
        X = [x for (x, _) in K]
        Y = [y for (_, y) in K]
        V = [v for (_, v) in band]
        return PlotlyJS.mesh3d(x=X, y=Y, z=V, intensity=V, colorscale="Viridis")
    end

    return PlotlyJS.plot([plotband(n|>getband) for n in 1:bandcount])
end

fiodir("/Users/alphaharrius/ZERData/squaremetal")
H = fioload("H")

energyspectrum = H|>crystalspectrum
groundstates = groundstatespectrum(energyspectrum, perunitcellfillings=0.5)

scale = Scale([4 0; 0 4], H|>getoutspace|>getcrystal|>getspace)
block = scale * getoutspace(H)
blockedH = block * H * block'

visualize(energyspectrum)

using LinearAlgebra

function Zipper.visualize(spectrum::CrystalSpectrum{1})
    eigenmodes = sort([(k|>norm, modes) for (k, modes) in spectrum|>geteigenmodes], by=(in -> in[1]))
    data = hcat(([spectrum.eigenvalues[mode] for mode in modes]|>sort for (_, modes) in eigenmodes)...)
    K = [v for (v, _) in eigenmodes]
    return plot(K, [data[n, :] for n in axes(data, 1)])
end

square = euclidean(RealSpace, 2)
point = [1/2, 1/2] ∈ square
spatialsnappingcalibration([point])

kspace = convert(MomentumSpace, square)

c4 = pointgrouptransform([0 -1; 1 0], localspace=square)
c2 = pointgrouptransform([-1 0; 0 -1], localspace=square)
m135 = pointgrouptransform([0 -1; -1 0], localspace=square)
m45 = pointgrouptransform([0 1; 1 0], localspace=square)

unitcell = Subset(point)
crystal = Crystal(unitcell, [128, 128])
reciprocalhashcalibration(crystal.sizes)

fiodir("/Users/alphaharrius/ZERData/squaremetal/RG1")
globaldistiller = fioload("rg1globaldistiller")
globaldistillerspectrum = globaldistiller|>crystalspectrum

frozenbands = distillation(globaldistillerspectrum, :frozen=>(v -> v > 0.5))[:frozen]
frozenbands|>visualize
frozenprojector = frozenbands|>crystalprojector
frozencorrelations = idmap(frozenprojector|>getoutspace) - frozenprojector

crystalfock = globaldistiller|>getoutspace
localcrystal = globaldistillerspectrum|>getcrystal
localspace = localcrystal|>getspace
normalvector = (1, 1) ∈ localspace

extendedrestrict = extendedcrystalrestrict(
        crystal=localcrystal, normalvector=normalvector, stripradius=0.5)
restrict = extendedrestrict * crystalfock
stripcorrelations = restrict * frozencorrelations * restrict'
stripspectrum = stripcorrelations|>crystalspectrum

stripspectrum|>linespectrum|>visualize

function Zipper.visualize(source::SnappingResult)
    p = scatter(
        LinRange(0, 1, source.denominator + 1), 
        zeros(source.denominator + 1), 
        markersize=3, markercolor=(:blue), markerstrokewidth=2, alpha=0.8, label="Snapping Grid")|>plot
    scatter!(
        source.forvalues, zeros(source.forvalues |> length), 
        markersize=7, markerstrokewidth=2, markercolor=(:orange), alpha=0.5, label="Input Values")
    return p
end

snappingdenominator(rand(5))|>visualize

function Zipper.visualize(source::EigenSpectrum)
    eigenfocks = FockSpace(m for (m, _) in source |> geteigenvalues) |> sparsegrouping(:eigenindex)
    datas::Vector = [[(source |> geteigenvalues)[m] |> real for m in subspace] for subspace in eigenfocks]
    sorted::Vector = [data for (_, data) in sort([data[1]=>data for data in datas], by=(v -> v.first))]
    paddinglengths::Vector = ([0, (data |> length for data in sorted)...]|>cumsum)[1:end-1]
    plotdatas = [
        ([paddinglength:(paddinglength + length(data) - 1)...], data)
        for (data, paddinglength) in zip(sorted, paddinglengths)]

    p = plot(scatter(plotdatas[1]..., markersize=5, markerstrokewidth=2, alpha=0.85), legend=false)
    for pdata in plotdatas[2:end]
       scatter!(pdata..., markersize=5, markerstrokewidth=2, alpha=0.85)
    end
    return p
end

r = blockedH.blocks|>first|>last|>eigspech|>visualize
r[1]
r[1]|>last
scatter(r[1]|>first, r[1]|>last)

""" Plot an `Element` for visualization. """
visualize(o::Element) = @warn "Plotting function for $(o|>typeof) is not defined."
export visualize

"""
    visualize(fockmap::FockMap)

Visualize the FockMap as a heatmap in real and imaginary part separately.
"""
function visualize(fockmap::FockMap)
    source = fockmap|>rep
    realplot = heatmap(source|>real, title="Re")
    imagplot = heatmap(source|>imag, title="Im")
    return plot(realplot, imagplot, yflip=true)
end

struct PlotRegion{D}
    region::Region
end

PlotRegion(region::Region) = PlotRegion{region|>getspace|>dimension}(region)

"""
    visualize(region::Region)

Visualize the region as a scatter plot.

### Keyword Arguments
- `visualspace` The `AffineSpace` to visualize the region.
"""
function visualize(region::Region; kwargs...)
    plotregion = PlotRegion(region)
    return visualize(plotregion, kwargs...)
end

""" The `2D` backend for plotting a region. """
function visualize(region::PlotRegion{2}; visualspace::AffineSpace = region.region|>getspace|>euclidean)
    r = [visualspace*r|>vec for r in region.region]
    x = [v|>first for v in r]
    y = [v|>last for v in r]
    # aspect_ratio=:equal for scatter is not working as expected, it only works when setting to 1.1.
    return scatter(x, y, aspect_ratio=1.1, xlims=(minimum(x)-1, maximum(x)+1), ylims=(minimum(y)-1, maximum(y)+1))
end

"""
    visualize(states::RegionState{2})

Visualize all the single particle state from a `RegionState` as a set of scatter plot, the size of each marker 
defines the magnitude of the function at that position; the color of the marker defines the complex phase of 
the function in `HSV` color scale. Hovering over each marker will show a truncated complex number of the function.
Each subplot of the state will have a title of `Mode {order}` and a table of the plotted modes will be shown in 
the terminal, and the table column `order` will corresponds to the title of the subplot.

### Keyword Arguments
- `visualspace` The `RealSpace` to visualize the region.
- `markersize`  The size of the marker.
- `logscale`    Depending on the function some marker will be much larger than the others, this parameter 
                will scale the size of the marker logarithmically, while a smaller value will decrease the 
                size distribution, a larger value will increase the size distribution. This value normally 
                should be between `0` and `1`.
"""
function visualize(
    states::RegionState{2};
    visualspace::RealSpace = states|>getoutspace|>getregion|>getspace|>euclidean,
    markersize::Integer = 1, logscale::Real = 1)

    function plotstate(order::Integer, statemap::SparseFockMap)
        r = [visualspace*getpos(mode)|>vec for mode in statemap|>getoutspace]
        x = [v|>first for v in r]
        y = [v|>last for v in r]
        f = ((statemap[mode, :]|>rep)[1, 1] for mode in statemap|>getoutspace)
        color = [convert(RGB{Float64}, HSV(angle(v) / 2π * 360, 1, 1)) for v in f]
        size = [max(markersize*(abs(v)^logscale), 1e-5) for v in f]
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

"""
    visualize(spectrum::CrystalSpectrum{1})

Visualize the eigenmodes of a `CrystalSpectrum` as a line plot, the y-axis is the spectrum values and 
the x-axis is sorted with the corresponding `Momentum` by their `norm`. Each band will be plotted as an 
individual line.
"""
function visualize(spectrum::CrystalSpectrum{1})
    eigenmodes = sort([(k|>norm, modes) for (k, modes) in spectrum|>geteigenmodes], by=(in -> in[1]))
    data = hcat(([spectrum.eigenvalues[mode] for mode in modes]|>sort for (_, modes) in eigenmodes)...)
    K = [v for (v, _) in eigenmodes]
    return plot(K, [data[n, :] for n in axes(data, 1)])
end

import PlotlyJS
"""
    visualize(spectrum::CrystalSpectrum{2})

Visualize the eigenmodes of a `CrystalSpectrum` as a 3D mesh plot, the z-axis is the spectrum values, 
the xy-plane will be the identity reciprocal space, and restricted to the brillouin zone of the crystal 
for plotting the spectrum. Each band will be plotted as an individual surface.

### Keyword Arguments
- `paddirection` The direction to pad the spectrum, either `:top` or `:bottom`, this is used when some bands 
                 at a `Momentum` have different count of eigenmodes than the others, the shorter bands will 
                 be padded with `NaN` to the direction specified, the `NaN` values will be ignored by `Plots.jl` 
                 and therefore shown as empty.  
"""
function visualize(spectrum::CrystalSpectrum{2}; paddirection::Symbol=:top)
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

"""
    visualize(source::SnappingResult)

Visualize the snapping result as a scatter plot, the blue markers specifies the snapping grid, and the orange
markers specifies the input values, from this plot you can observe how well the snapping is performing for your 
given values.
"""
function visualize(source::SnappingResult)
    p = scatter(
        LinRange(0, 1, source.denominator + 1), 
        zeros(source.denominator + 1), 
        markersize=3, markercolor=(:blue), markerstrokewidth=2, alpha=0.8, label="Snapping Grid")|>plot
    scatter!(
        source.forvalues, zeros(source.forvalues |> length), 
        markersize=7, markerstrokewidth=2, markercolor=(:orange), alpha=0.5, label="Input Values")
    return p
end

"""
    visualize(source::EigenSpectrum)

Visualize the eigenvalues of a `EigenSpectrum` as a scatter plot, the x-axis is the index of the eigen-modes, 
since the eigen-modes in the `EigenSpectrum` is grouped by their eigenvalues, each group of eigen-modes will 
have their own color in the plot.
"""
function visualize(source::EigenSpectrum)
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

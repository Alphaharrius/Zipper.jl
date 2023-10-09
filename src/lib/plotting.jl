function visualize(fockmap::FockMap; title::String = "", rowrange = :, colrange = :)
    source = rep(fockmap)[rowrange, colrange]
    layout = Layout(yaxis=attr(autorange="reversed"))
    subplots = [plot(heatmap(z=real(source)), layout) plot(heatmap(z=imag(source)), layout)]
    relayout!(subplots, title_text=title)
    subplots
end

function visualize(regions::Subset{<: Point}...; title::String = "", visualspace::AffineSpace = regions |> first |> getspace |> euclidean)
    function generateregionplot(region::Subset{<: Point})
        positions::Array{Vector} = [lineartransform(visualspace, point) |> vec for point in region]
        return visualizevectorpositions(positions)
    end
    # How to make all axis having the same scale:
    # https://stackoverflow.com/questions/52863305/plotly-scatter3d-how-can-i-force-3d-axes-to-have-the-same-scale-aspect-ratio
    layout = Layout(title=title, scene=attr(aspectmode="data"))
    fig = plot([region |> generateregionplot for region in regions], layout)
    relayout!(fig, template="simple_white")
    fig
end

function visualizevectorpositions(positions::Array{Vector})
    values::Matrix = hcat(positions...)
    @assert(0 < size(values, 1) <= 3)
    # Pad all positions to 3-dimensional to generalize all cases.
    padded_values::Matrix = zeros(3, size(values, 2))
    copyto!(view(padded_values, 1:size(values, 1), :), values)
    trace = scatter3d(x=padded_values[1,:], y=padded_values[2,:], z=padded_values[3,:], mode="markers")
    return trace
end

function visualize(spectrum::Vector{Pair{Mode, T}}; title::String = "") where {T <: Number}
    𝑁::Int64 = length(spectrum)
    ∑𝑝::Vector{Point} = [getattr(pair.first, :offset) + getattr(pair.first, :pos) for pair in spectrum]
    𝑀ₚ::Matrix{Float64} = hcat(map(𝑝 -> vec(lineartransform(euclidean(RealSpace, dimension(𝑝)), 𝑝)), ∑𝑝)...)
    markerpositions::Matrix{Float64} = zeros(3, 𝑁)
    copyto!(view(markerpositions, 1:size(𝑀ₚ, 1), :), 𝑀ₚ)
    sizes::Vector{Float64} = [abs(pair.second) for pair in spectrum]
    markersizes::Vector{Float64} = sizes / norm(sizes) * 120
    colors::Vector{RGB{Float64}} = [convert(RGB{Float64}, HSV(angle(pair.second) / 2π * 360, 1, 1)) for pair in spectrum]
    markercolors::Vector{Tuple{Float32, Float32, Float32, Float32}} = map(c -> Tuple([c.r, c.g, c.b, 1.0] * 255), colors)
    trace = scatter3d(
        x=markerpositions[1, :], y=markerpositions[2, :], z=markerpositions[3, :], mode="markers",
        marker=attr(
            symbol="circle",
            size=markersizes,
            color=markercolors))
    layout::Layout = Layout(title=title, scene=attr(aspectmode="data"))
    fig = plot([trace], layout)
    relayout!(fig, template="simple_white")
    fig
end

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

function visualize(spectrum::CrystalSpectrum{2}; title="", toppadding::Bool = true)
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
    plot([surface(z=plottingspectrum[n, :, :]) for n in axes(plottingspectrum, 1)], layout)
end

function visualize(source::SnappingResult; title::String = "")
    layout::Layout = Layout(title=title)
    plot([
        scatter(x=LinRange(0, 1, source.denominator + 1), y=zeros(source.denominator + 1), mode="markers", name="markers", marker=attr(
            symbol=:square,
            line=attr(width=2, color="midnightblue"),
            color="lightskyblue", size=10
        )),
        scatter(x=source.forvalues, y=zeros(source.forvalues |> length), mode="markers", name="values")], layout)
end

function visualize(source::EigenSpectrum; title::String = "")
    eigenfock::FockSpace = FockSpace(m for (m, _) in source |> geteigenvalues) |> sparsegrouping(:eigenindex)
    realdatas::Vector = [[(source |> geteigenvalues)[m] |> real for m in subspace] for subspace in eigenfock |> subspaces]
    # TODO: Think of a better way to visualize complex datas.

    function generateplot(datas::Vector, symbol::Symbol)::Vector
        sorteddatas::Vector = [data for (_, data) in sort([data[1] => data for data in datas], by=v -> v.first)]
        paddinglengths::Vector = ([0, (data |> length for data in sorteddatas)...] |> cumsum)[1:end-1]
        return [scatter(x=paddinglength:(paddinglength + length(data)), y=data, mode="markers", marker=attr(size=10, line_width=2, symbol=symbol)) for (data, paddinglength) in zip(sorteddatas, paddinglengths)]
    end

    layout = Layout(title=title)
    plot([generateplot(realdatas, :circle)...], layout)
end

export visualize

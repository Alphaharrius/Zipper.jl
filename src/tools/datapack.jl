mutable struct RGData
    data::Dict{Symbol, Any}
    headloc::Vector{String}
    head::Dict{Symbol, Any}
end
export RGData

function Base.:show(io::IO, v::RGData)
    headlocname = join(v.headloc, "->")
    print(io, "$(v|>typeof)($headlocname)")
end

function rgpack(correlations)
    data = Dict{Symbol, Any}(:target=>correlations)
    return RGData(data, [], data)
end
export rgpack

function getstepper(;f::Function, name::String)
    function stepper!(data::RGData)
        target = data.head[:target]

        kws = Base.kwarg_decl(f|>methods|>first)
        !isempty(kws) && !haskey(data.head, :computeds) && error("No computed data available...")
        kwargs = (kw=>data.head[:computeds][kw|>String] for kw in kws)

        stepped, step = isempty(kws) ? f(target) : f(target; kwargs...)
        steps = haskey(data.head, :steps) ? data.head[:steps] : (data.head[:steps] = Dict{String, Any}())
        steps[name] = Dict{Symbol, Any}(:target=>stepped, :step=>step)
        data.head = steps[name]
        push!(data.headloc, name)
        return data
    end
    return stepper!
end
export getstepper

function getcomputator(;f::Function, name::String)
    function computator!(data::RGData)
        target = data.head[:target]

        kws = Base.kwarg_decl(f|>methods|>first)
        !isempty(kws) && !haskey(data.head, :computeds) && error("No computed data available...")
        kwargs = (kw=>data.head[:computeds][kw|>String] for kw in kws)

        args = Base.method_argnames(f|>methods|>first)

        computed = isempty(kws) ? f(target) : length(args) > 1 ? f(target; kwargs...) : f(; kwargs...)
        computeds = haskey(data.head, :computeds) ? data.head[:computeds] : (data.head[:computeds] = Dict{String, Any}())
        computeds[name] = computed
        return data
    end
    return computator!
end
export getcomputator

gettarget(data::RGData) = data.head[:target]
export gettarget

getcomputed(data::RGData, name) = data.head[:computeds][name|>String]
export getcomputed

function tracedata(d::Dict, subkeys)
    target = d
    for key in subkeys
        target = target[:steps][key]
    end
    return target
end

function getsteps(data::RGData)
    if !haskey(data.head, :steps)
        return []
    end
    return data.head[:steps]|>keys|>collect
end
export getsteps

function revert!(data::RGData)
    if isempty(data.headloc)
        error("Cannot revert further...")
    end
    pop!(data.headloc)
    data.head = tracedata(data.data, data.headloc)
    return data
end
export revert!

function step!(data::RGData, name::String)
    name in getsteps(data) || error("Step $name not found...")
    push!(data.headloc, name)
    data.head = tracedata(data.data, data.headloc)
    return data
end
export step!

function head!(data::RGData, names...)
    head = tracedata(data.data, names)
    data.headloc = names|>collect
    data.head = head
    return data
end
export head!

function base!(data::RGData)
    data.head = data.data
    data.headloc = []
    return data
end
export base!

function getsteptransform(data::RGData, names...)
    transform = 1
    current = data.head
    for name in names
        haskey(current, :steps) || error("Cannot go any further to $name...")
        haskey(current[:steps], name) || error("Step $name not found...")
        transform = transform*current[:steps][name][:step]
        current = current[:steps][name]
    end
    return transform
end
export getsteptransform

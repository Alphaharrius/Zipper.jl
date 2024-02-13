MEMORIZED::Dict{Symbol, Dict} = Dict()

function memorization(f::Function)
    fsymbol = Symbol(f)
    haskey(MEMORIZED, fsymbol) || (MEMORIZED[fsymbol] = Dict())
    return function (arguments...; kwarguments...)
        key = Dict(:arguments=>arguments, kwarguments...)
        if haskey(MEMORIZED[fsymbol], key)
            @debug "Using memorized value for $(f) with arguments $(arguments)"
            return MEMORIZED[fsymbol][key]
        else
            result = f(arguments...; kwarguments...)
            MEMORIZED[fsymbol][key] = result
            return result
        end
    end
end

macro memorize(innerdef)
    innerdef = ExprTools.splitdef(innerdef)
    outerdef = copy(innerdef)
    fname = get(innerdef, :name, nothing)
    if fname !== nothing
        @assert fname isa Symbol
        innerdef[:name] = Symbol(fname, :_unmemorized)
    end
    outerdef[:body] = Expr(:call,
        :($memorization($(ExprTools.combinedef(innerdef)))),
        get(outerdef, :args, [])...,
        get(outerdef, :kwargs, [])...,
    )
    return esc(ExprTools.combinedef(outerdef))
end
export @memorize

MEMORIZED::Dict{Symbol, Dict} = Dict()

function memoization(f::Function)
    fsymbol = Symbol(f)
    haskey(MEMORIZED, fsymbol) || (MEMORIZED[fsymbol] = Dict())
    return function (arguments...; kwarguments...)
        key = Dict(:arguments=>arguments, kwarguments...)
        if haskey(MEMORIZED[fsymbol], key)
            @debug "Using memoized value for $(f) with arguments $(arguments)"
            return MEMORIZED[fsymbol][key]
        else
            result = f(arguments...; kwarguments...)
            MEMORIZED[fsymbol][key] = result
            return result
        end
    end
end

macro memoize(innerdef)
    innerdef = ExprTools.splitdef(innerdef)
    outerdef = copy(innerdef)
    fname = get(innerdef, :name, nothing)
    if fname !== nothing
        @assert fname isa Symbol
        innerdef[:name] = Symbol(fname, :_default)
    end
    outerdef[:body] = Expr(:call,
        :($memoization($(ExprTools.combinedef(innerdef)))),
        get(outerdef, :args, [])...,
        get(outerdef, :kwargs, [])...,
    )
    return esc(ExprTools.combinedef(outerdef))
end
export @memoize

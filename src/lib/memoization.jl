MEMOIZED::ConcurrentDict{Symbol, ConcurrentDict} = ConcurrentDict{Symbol, ConcurrentDict}()

function memoization(f::Function)
    fsymbol = Symbol(f)
    haskey(MEMOIZED, fsymbol) || (MEMOIZED[fsymbol] = ConcurrentDict{Dict, Any}())
    return function (arguments...; kwarguments...)
        key = Dict(:arguments=>arguments, kwarguments...)
        if haskey(MEMOIZED[fsymbol], key)
            @debug "Using memoized value for $(f) with arguments $(arguments)"
            return MEMOIZED[fsymbol][key]
        else
            result = f(arguments...; kwarguments...)
            MEMOIZED[fsymbol][key] = result
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

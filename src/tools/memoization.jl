MEMOIZED::ConcurrentDict{Symbol, ConcurrentDict} = ConcurrentDict{Symbol, ConcurrentDict}()

abstract type AbstractArguments end

struct Arguments <: AbstractArguments
    arguments::Tuple
end

struct FullArguments <: AbstractArguments
    arguments::Tuple
    kwarguments::Dict
end

Base.:hash(o::Arguments)::UInt = reduce(⊻, argument|>hash for argument in o.arguments)

Base.:hash(o::FullArguments)::UInt = reduce(⊻, argument|>hash for argument in o.arguments) ⊻ hash(o.kwarguments)

Base.:(==)(a::Arguments, b::Arguments)::Bool = a.arguments == b.arguments

Base.:(==)(a::FullArguments, b::FullArguments)::Bool = a.arguments == b.arguments && a.kwarguments == b.kwarguments

function memoization(f::Function)
    fsymbol = Symbol(f)
    haskey(MEMOIZED, fsymbol) || (MEMOIZED[fsymbol] = ConcurrentDict{AbstractArguments, Any}())
    return function (arguments...; kwarguments...)
        key = length(kwarguments) == 0 ? Arguments(arguments) : FullArguments(arguments, Dict(kwarguments))
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

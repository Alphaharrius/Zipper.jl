MEMORIZED::Dict{Symbol, Dict} = Dict()

function memorization(f::Function)
    fsymbol = Symbol(f)
    haskey(MEMORIZED, fsymbol) || (MEMORIZED[fsymbol] = Dict())
    return function (arguments...)
        key = Tuple(arguments...)
        if haskey(MEMORIZED[fsymbol], key)
            @info "Using memorized value for $(f) with arguments $(key)"
            return MEMORIZED[fsymbol][key]
        else
            result = f(arguments...)
            MEMORIZED[fsymbol][key] = result
            return result
        end
    end
end

macro memorize(f)
    :(memorization($f))
end
export @memorize

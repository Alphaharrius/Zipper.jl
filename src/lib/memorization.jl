MEMORIZED::Dict{Symbol, Dict} = Dict()

function memorization(f::Function)
    fsymbol = Symbol(f)
    haskey(MEMORIZED, fsymbol) || (MEMORIZED[fsymbol] = Dict())
    return function (arguments...)
        if haskey(MEMORIZED[fsymbol], arguments)
            @info "Using memorized value for $(f) with arguments $(arguments)"
            return MEMORIZED[fsymbol][arguments]
        else
            result = f(arguments...)
            MEMORIZED[fsymbol][arguments] = result
            return result
        end
    end
end

macro memorize(f)
    :(memorization($f))
end
export @memorize

# LZW compression implementation from https://rosettacode.org/wiki/LZW_compression#Julia

function lzwcompress(decompressed::String)::Vector{Int32}
    dictsize = 256
    dict     = Dict{String,Int}(string(Char(i)) => i for i in 0:dictsize)
    result   = Vector{Int32}(undef, 0)
    w        = ""
    watchprogress(desc="lzwcompress")
    for c in decompressed
        wc = string(w, c)
        if haskey(dict, wc)
            w = wc
        else
            push!(result, dict[w])
            dict[wc]  = dictsize
            dictsize += 1
            w        = string(c)
        end
        updateprogress()
    end
    unwatchprogress()
    if !isempty(w) push!(result, dict[w]) end
    return result
end
 
function lzwdecompress(compressed::Vector{Int32})
    dictsize = 256
    dict     = Dict{Int,String}(i => string('\0' + i) for i in 0:dictsize)
    result   = IOBuffer()
    w        = string(Char(popfirst!(compressed)))
    write(result, w)
    watchprogress(desc="lzwdecompress")
    for k in compressed
        if haskey(dict, k)
            entry = dict[k]
        elseif k == dictsize
            entry = string(w, w[1])
        else
            error("bad compressed k: $k")
        end
        write(result, entry)
        dict[dictsize] = string(w, entry[1])
        dictsize += 1
        w = entry
        updateprogress()
    end
    unwatchprogress()
    return String(take!(result))
end

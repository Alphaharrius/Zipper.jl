using Zipper

# ▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃
# ◆  codeformatting.jl ◆
# This file is not supposed to be included to the base library since this is an 
# embedded formatter for the source code which will add some decorations to make 
# it more readable and understandable. This formatter will recognise the prefix 
# in the comments to inject pretty splitters and titles to the source code without 
# affecting the functionality of the code.
# ▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃

# ▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃
# ◆  Modified definition ◆
struct Modified
    path::String
    lines::Vector{String}
    updatedcount::Int
end
# ▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃

# ▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃
# ◆ Modified display ◆
Base.:show(io::IO, m::Modified) = print(io, "$(m|>typeof)(updated=$(m.updatedcount))")
# ▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃

# ▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃
# ◆ Modified APIs ◆
function formatsource(path::String)
    lines::Vector{String} = []
    maxlinelength::Int = 0
    open(path, "r") do sourcefile
        for line in eachline(sourcefile)
            maxlinelength = max(maxlinelength, line|>length)
            push!(lines, line)
        end
    end
    beginsyntax = "#\\begin:"
    endsyntax = "#\\end"
    splitterline = "# "*"▃"^(maxlinelength-2)
    Zipper.watchprogress(desc="formatsource")
    updatedcount::Int = lines|>length
    processed = []
    for line in lines
        if startswith(line, beginsyntax)
            push!(processed, splitterline)
            push!(processed, "# ◆ "*line[length(beginsyntax)+1:end]*" ◆")
        elseif startswith(line, endsyntax)
            push!(processed, splitterline)
        else
            push!(processed, line)
            updatedcount -= 1
        end
        Zipper.updateprogress()
    end
    Zipper.unwatchprogress()
    return Modified(path, processed, updatedcount)
end

function checklines(modified::Modified)
    for line in modified.lines
        println(line)
    end
end

function applymodifications(modified::Modified)
    open(modified.path, "w") do sourcefile
        for line in modified.lines
            println(sourcefile, line)
        end
    end
end
# ▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃

modified = formatsource(
    "/Users/alphaharrius/Code Projects/Zipper.jl/src/lib/transformations.jl")
checklines(modified)
applymodifications(modified)

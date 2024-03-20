#\begin:Collection tools
"""
    matchreduce(left, right; leftmatcher::Function, rightmatcher::Function, reducer::Function)

Matches elements from two collections `left` and `right` based on provided matching functions, 
and then applies a reducing function to each pair of matched elements.

### Input
- `left`        The first collection.
- `right`       The second collection.
- `leftmatch`   Function to generate keys for matching elements from `left`. This function 
                should take an element from `left` and return a key that will be used to find 
                a matching element in `right`.
- `rightmatch`  Function to generate keys for matching elements from `right`. This function 
                should take an element from `right` and return a key that will be used to create 
                a lookup dictionary.
- `reducer`     Function to apply to each pair of matched elements, with two arguments which the 
                first will take in the matched element of the `left` collection; and the second 
                will take in the matched element of the `right` collection.

### Output
A generator of reduced values for each pair of matched elements.
"""
function matchreduce(
    left, right; leftmatch::Function, rightmatch::Function, reducer::Function = (a, b)->(a, b))

    lookup::Dict = Dict((v|>rightmatch)=>v for v in right)
    leftmatchings = ((v, v|>leftmatch) for v in left)
    valid = ((v, k) for (v, k) in leftmatchings if haskey(lookup, k))
    matched::Base.Generator = ((v, lookup[k]) for (v, k) in valid)

    return (reducer(l, r) for (l, r) in matched)
end

function bijection(o::Dict)::Dict
    bijected = Dict(v => k for (k, v) in o)
    length(o) == length(bijected) || @error "The mappings are not bijective!"
    return bijected
end
#\end

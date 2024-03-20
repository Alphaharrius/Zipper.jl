# ▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃
# ◆ Mode definition ◆
"""
    Mode(attrs::Dict{Symbol})
    Mode(generator::Base.Generator)
    Mode(input::Base.Iterators.Flatten)
    Mode(datas::Vector{Pair{Symbol, T}}) where {T}

Represents an element in a `FockSpace`, and uniquely identifies a physical mode. 
The attribute `:orbital` is defaulted to `swave` if not specified.

### Attributes to put in `attrs`
- `:r` stores a `Offset` which is the offset in lattice unit.
- `:k` stores a `Momentum` which is the offset in reciprocal lattice unit.
- `:b` stores a `Point` which is the unit cell offset, this is associated to the attribute `:flavor`.
- `:flavor` stores an `Integer` that identifies a fermionic freedom at a lattice site.
- `:orbital` Defines which orbital this mode transforms like under a symmetry, for example for 
             𝐶₃ symmetry and 𝑠 like orbital `Dict(:c3 => :s)`.

### Input
- `attrs` The attributes which uniquely identifies the `Mode` object.
"""
struct Mode <: Element{Mode}
    attrs::Dict{Symbol}

    Mode(attrs::Dict{Symbol}) = new(Dict((:orbital => swave, attrs...)))

    Mode(generator::Base.Generator) = Mode(Dict(generator))
    Mode(input::Base.Iterators.Flatten) = Mode(Dict(input))

    Mode(attrs::Vector{Pair{Symbol, T}}) where {T} = Mode(Dict(attrs...))
    Mode(attrs::Pair{Symbol, <:Any}...) = Mode(Dict(attrs...))
end
export Mode
# ▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃

# ▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃
# ◆ Mode essentials ◆

# Mode can be used as keys in dictionaries and sets.
Base.:hash(mode::Mode)::UInt = hash(mode.attrs)
# Since Mode is a wrapper around a Dict, supporting getting all attribute keys.
Base.keys(mode::Mode) = mode.attrs|>keys
# ▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃

# ▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃
# ◆ Mode logicals ◆
Base.:(==)(a::Mode, b::Mode)::Bool = a.attrs == b.attrs
# ▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃

# ▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃
# ◆ Mode arithmetics ◆
Base.:+(mode::Mode, offset::Offset)::Mode = setattr(mode, :r => getattr(mode, :r) + offset)
Base.:+(mode::Mode, offset::Momentum)::Mode = setattr(mode, :k => getattr(mode, :k) + offset)
Base.:-(mode::Mode, offset::Point)::Mode = mode + (-offset)
# ▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃

# ▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃
# ◆ Mode attribute API ◆
"""
    hasattr(mode::Mode, key::Symbol)::Bool

Check if the `mode` has the attribute identified by `key`.
"""
hasattr(mode::Mode, key::Symbol)::Bool = haskey(mode.attrs, key)
export hasattr

"""
    getattr(mode::Mode, key::Symbol)

Retrieve the attribute value identified by `key` from `mode`.
"""
getattr(mode::Mode, key::Symbol) = mode.attrs[key]
export getattr

"""
    getattr(key::Symbol)

Shorthand of `getattr(v, key::Symbol)` with the pipe operator `|>`.

### Examples
The line `mode |> getattr(:flavor)` is equal to `getattr(mode, :flavor)`.
"""
getattr(key::Symbol) = v -> getattr(v, key)

""" A shorthand for getting the mode attributes for visualization in the REPL. """
getattrs(mode::Mode) = mode.attrs
export getattrs

"""
    removeattr(mode::Mode, keys::Symbol...)::Mode

Create a **copy** of `mode` **without** the attributes identified by `keys`.

### Examples
- To remove the attribute of `:r` and `:b`, we use `removeattr(mode, :r, :b)`.
"""
removeattr(mode::Mode, keys::Symbol...)::Mode = Mode(Dict(filter(p -> !(p.first ∈ keys), mode.attrs)))
export removeattr

"""
    removeattr(keys::Symbol...)

Shorthand of `removeattr(v, keys::Symbol...)::Mode` with the pipe operator `|>`.

### Examples
The line `mode |> removeattr(:flavor, :index)` is equal to `removeattr(mode, :flavor, :index)`.
"""
removeattr(keys::Symbol...) = v -> removeattr(v, keys...)

"""
    removeattr(modes, keys::Symbol...)::Subset{Mode}

Create a **copy** of every `Mode` of `modes` **without** the attributes identified by `keys`, the 
resulting `Subset` might not have the same length as the input `modes` as some `Mode` might be **condensed** 
into a single one after some unique identifier attributes is removed.

### Examples
To remove the attribute of `:r` and `:b`, we use `removeattr(modes, :r, :b)`.
"""
removeattr(modes, keys::Symbol...)::Subset{Mode} = Subset(removeattr(mode, keys...) for mode in modes)

"""
    setattr(mode::Mode, attrs::Pair{Symbol}...)::Mode

Create a **copy** of `mode` with the new attributes identified by `keys` added to the current attributes, 
if the attribute exists in `mode`, the current record will be overwritten.

### Examples
- To add `:r` and `:flavor`, we use `newmode = setattr(mode, :r => point, :flavor => 1)`
"""
setattr(mode::Mode, attrs::Pair{Symbol}...)::Mode = Mode(Dict(mode.attrs..., attrs...))
export setattr

"""
    setattr(attrs::Pair{Symbol}...)

Shorthand of `setattr(v, attrs::Pair{Symbol}...)::Mode` with the pipe operator `|>`.

### Examples
The line `mode |> setattr(:flavor => 1)` is equal to `setattr(mode, :flavor => 1)`.
"""
setattr(attrs::Pair{Symbol}...) = v -> setattr(v, attrs...)

""" Allow merging attributes of two distinct `Mode`. """
Base.:merge(a::Mode, b::Mode)::Mode = Mode(Dict(getattrs(a)..., getattrs(b)...))

""" A shorthand for setting an attribute to a mode. """
Base.:&(mode::Mode, attr::Pair{Symbol}) = setattr(mode, attr)

"""
    setattr(subset, attrs::Pair{Symbol}...)::Subset{Mode}

Create a **copy** of `modes` with the new attributes identified by `keys` added to the current attributes 
for each mode in `modes`, if the attribute exists in the modes, the current record will be overwritten.

### Examples
- To add `:r` and `:flavor`, we use `newmodes = setattr(modes, :r => point, :flavor => 1)`
"""
setattr(subset, attrs::Pair{Symbol}...)::Subset{Mode} = Subset(setattr(mode, attrs...) for mode in subset)

"""
    renameattr(mode::Mode, oldkey::Symbol, newkey::Symbol)::Mode

Renames an attribute of a `Mode` object.

# Arguments
- `mode` The `Mode` object whose attribute is to be renamed.
- `oldkey` The current name of the attribute.
- `newkey` The new name for the attribute.

# Returns
- `Mode`: A new `Mode` object with the renamed attribute.
"""
renameattr(mode::Mode, oldkey::Symbol, newkey::Symbol)::Mode = (
    mode|>setattr(newkey=>getattr(mode, oldkey))|>removeattr(oldkey))
renameattr(input, oldkey::Symbol, newkey::Symbol) = Subset(
    renameattr(mode, oldkey, newkey) for mode in input)
renameattr(oldkey, newkey)::Function = m -> renameattr(m, oldkey, newkey)

function renameattr(renamings::Pair{Symbol, Symbol}...)::Function
    function rename(mode::Mode)::Mode
        ret::Mode = mode
        for (oldkey, newkey) in renamings
            ret = renameattr(ret, oldkey, newkey)
        end
        return ret
    end
    return rename
end

export renameattr

"""
    commonattr(modes, key::Symbol)

Retrieve the common attribute associated to `key` of all the child modes of the `fockspace`, 
and throws assertion error if the attribute is not unique within the `fockspace`.
"""
function commonattr(modes, key::Symbol)
    set::Set = Set()
    foreach(m -> push!(set, getattr(m, key)), modes)
    @assert(length(set) == 1, "The modes in this fockspace does not share the same attr `$(key)`!")
    return first(set)
end
export commonattr

"""
    commonattrs(modes)::Base.Generator

Retrieve all the attributes within the group of modes that has the same values and returned 
as a generator of `Symbol`.
"""
function commonattrs(modes)::Base.Generator
    attrs::Base.Generator = (Set(k => v for (k, v) in mode |> getattrs) for mode in modes)
    return (k for (k, _) in intersect(attrs...))
end
export commonattrs

""" Get all the available attribute keys of all the modes. """
attrkeys(modes) = Set(key for mode in modes for key in mode|>keys)

"""
    indexmodes(modes)::Subset{Mode}

Retrieve the representative modes of this group of modes with their 
corresponding common attributes stripped.
"""
function indexmodes(modes)::Subset{Mode}
    cleanattrs::Base.Generator = modes |> commonattrs
    return Subset(mode for mode in modes) |> removeattr(cleanattrs...) 
end
export indexmodes

getorbital(mode::Mode, default::BasisFunction)::BasisFunction = (
    hasattr(mode, :orbital) ? getattr(mode, :orbital) : default)
export getorbital

getorbital(default::BasisFunction = swave) = mode -> getorbital(mode, default)

setorbital(mode::Mode, basis::BasisFunction)::Mode = setattr(mode, :orbital => basis)
export setorbital

setorbital(basis::BasisFunction) = mode -> setorbital(mode, basis)

"""
    mapmodes(mapper::Function, modes)::Subset{Mode}

A standard way to apply a mapping function onto an iterator of `Mode` that might 
ends up with mode duplication. All duplications will be mapped with incremental 
`:flavor` indices starting from `1`. It is adviced that all code that looks like 
`Subset(mode|>mapper for mode in modes)` to be replaced with this function.

### Input
- `mapper`  The mapping function to be applied to each `Mode` object, it should 
            only be taking one `Mode` as argument.
- `modes`   The iterator of `Mode` objects to be mapped.
"""
function mapmodes(mapper::Function, modes)::Subset{Mode}
    # :flavor is removed to ensure there is not overlooked degeneracy lifted.
    mapped::Base.Generator = (mode|>mapper|>removeattr(:flavor) for mode in modes)
    degenerates::Dict = Dict()
    flavoredmodes::Vector{Mode} = []
    for mode in mapped
        if !haskey(degenerates, mode)
            degenerates[mode] = 1
        else
            degenerates[mode] += 1
        end
        push!(flavoredmodes, mode|>setattr(:flavor=>degenerates[mode]))
    end
    return Subset(flavoredmodes)
end
export mapmodes

""" A shorthand for `mapmodes`. """
mapmodes(mapper::Function)::Function = input -> mapmodes(mapper, input)
# ▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃

# ▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃
# ◆ Mode display ◆
""" Display the number of attributes that identifies this `Mode`. """
Base.:show(io::IO, mode::Mode) = print(io, string("$(typeof(mode))$(tuple(keys(mode.attrs)...))"))

""" Shorthand for accessing the attribute information of all modes within the `FockSpace`. """
modeattrs(modes)::OrderedSet{Dict} = OrderedSet(mode|>getattrs for mode in modes)
export modeattrs

"""
    showmodes(modes; showkeys)

Display the mode attribute information in form of a table.

### Input
- `modes`    The modes to be displayed.
- `showkeys` The keys to be displayed, if not provided, all keys will be displayed.
"""
function showmodes(modes; showkeys=modes|>attrkeys|>collect)
    header = ["order", showkeys...]
    getrow(mode) = (hasattr(mode, key) ? mode|>getattr(key) : "N/A" for key in showkeys)
    rowlength = length(header)
    rows = (reshape([n, getrow(mode)...], (1, rowlength)) for (n, mode) in modes|>enumerate)
    pretty_table(vcat(rows...), header=header)
end
export showmodes
# ▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃

# ▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃
# ◆ Spatial API implementation ◆
"""
    getspace(mode::Mode)

The space of a `Mode` comes from the physical quantities its defined on, such as 
`:r`, `:k` and `:b`, if none of those are defined, it will be `euclidean(RealSpace, 1)` 
as the mode and it's siblings can always be parameterized by a scalar.

### Output
The space of the attribute `:k` or `:r` if available, fall back to `:b` if otherwise; 
returns a Euclidean space of dimension `1` no position attributes is found.
"""
function Zipper.:getspace(mode::Mode)
    # :r and :k have a higher priority in determining the space of the mode.
    if hasattr(mode, :k) return getspace(getattr(mode, :k)) end
    if hasattr(mode, :r) return getspace(getattr(mode, :r)) end
    if hasattr(mode, :b) return getspace(getattr(mode, :b)) end
    # If the mode does not based on any physical position or quantities for 
    # associating with any space, then it will simply lives in a 1D euclidean 
    # space as the mode and it's siblings can always be parameterized by a scalar.
    return euclidean(RealSpace, 1)
end

"""
    getpos(mode::Mode)::Point

Get the actual position of the mode, if the mode is associated with a `Momentum` which a `:k` 
attribute is attached, the return value will be the `Momentum`; if the mode is associated with 
a `Offset` which a `:r` attribute is attached, it will be returned, and if there is a intra 
unit cell offset `:b` attached, it will be added to the `Offset` as the final result.
"""
function Zipper.:getpos(mode::Mode)::Point
    if hasattr(mode, :k)
        return getattr(mode, :k)
    end

    if hasattr(mode, :r)
        return hasattr(mode, :b) ? getattr(mode, :r) + getattr(mode, :b) : getattr(mode, :r)
    end

    error("The mode does not have a primary position attribute!")
end
# ▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃

# ▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃
# ◆ API extensions ◆
""" Added to support shorthand of `Subset(m0, m1, ...)`. """
Zipper.:Subset(modes::Mode...) = Subset(modes)

"""
    spanoffset(basismodes::Subset{Mode}, points::Subset{<: Point})::Subset{Mode}

Given a set of `basismodes`, and the generator `points`, span the basis modes with 
the primary position attribute `:r` or `:k`, the primary ordering will be the ordering 
of `points`, then follows the ordering of `basismodes`.
"""
spanoffset(basismodes::Subset{Mode}, momentums::Subset{Momentum})::Subset{Mode} = Subset(
    setattr(mode, :k=>k) for k in momentums for mode in basismodes)
spanoffset(basismodes::Subset{Mode}, offsets::Subset{Offset})::Subset{Mode} = Subset(
    setattr(mode, :r=>point) for point in offsets for mode in basismodes)
export spanoffset
# ▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃

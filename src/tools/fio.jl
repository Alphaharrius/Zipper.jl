"""
### Attributes
- `projectpath` The path to the project directory which all object save/load will be performed.
- `threadtargetname` The name of the current FIO target of each thread in the Julia environment, 
                     this have took advantage of each thread can only save one object at a time.
"""
mutable struct FioState
    projectpath::String
    threadtargetname::Vector{String}
end

global FIO_STATE = FioState(
    pwd(), # The default directiory is the current working directory.
    ["" for i in 1:Threads.nthreads()]) # Defaults to an empty string.

"""
    fiodir(path::String)

Set the project directory to the given `path`, which all object save/load will be performed. 
If the directory does not exist, the function will attempt to create it. If the directory cannot be created, 
an error message will be logged.

### Input
- `path::String`: The path to the directory to be set as the resource directory.
"""
function fiodir(path::String)
    if !isdir(path)
        try
            mkdir(path)
        catch _
            error("Could not create resource directory at $path")
        end
    end
    FIO_STATE.projectpath = path
    @info "Resource directory is set to $path"
    return
end
export fiodir

"""
    fiodir()

Return the current project directory path.
"""
fiodir() = FIO_STATE.projectpath

# The JSON key of the type of the object.
const JL_TYPE::String = ":t"
# The JSON key of the data of the object.
const JL_DATA::String = ":d"

# A Dict that stores the special lower functions for the serialization of the object attributes.
global SPECIAL_LOWERS::Dict = Dict()

"""
    fiolower(typeattr::Tuple{Type, Symbol}, lower::Function)

Define a special lower function for the serialization of the object attributes, the lower function will 
be applied to the object before `JSON` serializes the object into a JSON string.

### Input
- `typeattr` The attribute of the object to be serialized, in format of `(::Type, ::Symbol)` 
             where the `::Symbol` is the name of the attribute.
- `lower` The lower function to be applied to the attribute value before serialization.
"""
fiolower(typeattr::Tuple{Type, Symbol}, lower::Function) = SPECIAL_LOWERS[typeattr] = lower
export fiolower

# A Dict that stores the special constructor functions for the deserialization of the objects.
global SPECIAL_CONSTRUCTORS::Dict = Dict()

"""
    fioconstructor(type::Type, constructor::Function)

Define a special constructor function for the deserialization of the object, the constructor function will 
be called after all object attributes are deserialized.

### Input
- `type` The type of the object to be deserialized.
- `constructor` The constructor function, noted that the constructor must have parameters in the order of 
                the actual ordering of attributes of the `type`.
"""
fioconstructor(type::Type, constructor::Function) = SPECIAL_CONSTRUCTORS[type] = constructor
export fioconstructor

# A Dict that stores the special parser functions for the deserialization of the object attributes.
global SPECIAL_PARSERS::Dict = Dict()

"""
    fioparser(typeattr::Tuple{Type, Symbol}, parser::Function)

Define a special parser function for the deserialization of the object attributes, the parser function will
be applied to the attribute value after the object is deserialized into a `Dict`.

### Input
- `typeattr` The attribute of the object to be deserialized, in format of `(::Type, ::Symbol)` 
             where the `::Symbol` is the name of the attribute.
- `parser` The parser function to be applied to the attribute value after deserialization.
"""
fioparser(typeattr::Tuple{Type, Symbol}, parser::Function) = SPECIAL_PARSERS[typeattr] = parser
export fioparser

"""
A wrapper to the `JSON` library for a custom lowering method for the `<:Element` types, 
the output will contain the Julia type information that assist in later deserializations.
"""
function JSON.lower(o::Union{Element, RGData})
    type::Type = typeof(o)
    data = Dict()
    for fieldname in type|>fieldnames
        fieldvalue = getfield(o, fieldname)
        # Parametric types will have the same structure of its original wrapper type.
        key = (type.name.wrapper, fieldname)
        if haskey(SPECIAL_LOWERS, key)
            data[fieldname] = SPECIAL_LOWERS[key](fieldvalue)
        else
            data[fieldname] = fieldvalue
        end
        updateprogress()
    end
    return Dict(JL_TYPE=>type|>string, JL_DATA=>data)
end

"""
Deserialize the JSON Dict of Zipper types back into the original Julia objects.
"""
function fiozipperparse(type::Type, data::Dict)
    arguments = []
    for fieldname in type|>fieldnames
        key = (type.name.wrapper, fieldname)
        # Parametric types will have the same structure of its original wrapper type.
        parser::Function = haskey(SPECIAL_PARSERS, key) ? SPECIAL_PARSERS[key] : fioparse
        fieldvalue = data[fieldname|>string]|>parser
        push!(arguments, fieldvalue)
        updateprogress()
    end
    constructor = type.name.wrapper
    haskey(SPECIAL_CONSTRUCTORS, type.name.wrapper) && (constructor = SPECIAL_CONSTRUCTORS[type.name.wrapper])
    haskey(SPECIAL_CONSTRUCTORS, type) && (constructor = SPECIAL_CONSTRUCTORS[type])
    return constructor(arguments...)
end

# General types will be passed.
fioparse(o::Any) = o

"""
Identity and parse the serialized Zipper types recursively from the JSON Dict deserialized by `JSON.parse`.
"""
function fioparse(object::Dict)
    if haskey(object, JL_TYPE)
        objecttype::Type = Meta.parse(object[JL_TYPE])|>eval
        return fiozipperparse(objecttype, object[JL_DATA])
    end
    processed::Dict = Dict()
    for (k, v) in object
        processed[k] = fioparse(v)
        updateprogress()
    end
    return processed
end

# Vectors will have their elements being parsed.
fioparse(vector::Vector) = [fioparse(el) for el in vector]

# A table to record the storage type of the specific Julia type in format of type=>storagetype.
global STORAGE_TYPES::Dict = Dict()
# A table to record the parsed type of the specific storage type in format of storagetype=>type.
global PARSED_TYPES::Dict = Dict()

"""
    fiostoragetype(type::Type, storage::Type)

Register a storage type for the given Julia type, the storage type will be serialized instead of the actual 
type when saving the object to a file. During deserialization the storage type will be deserialized and converted 
back to the actual Julia type. Please be noted that the storage type must be a subtype of `Element` and have 
the function `convert(::Type{StorageType}, ::JuliaType)` and `convert(::Type{JuliaType}, ::StorageType)` defined.

### Input
- `type` The Julia type to be registered.
- `storage` The storage type to be registered.
"""
function fiostoragetype(type::Type, storage::Type)
    STORAGE_TYPES[type] = storage
    PARSED_TYPES[storage] = type
    @debug "Added storage type $storage for $type"
end
export fiostoragetype

"""
    fiotargetname()

Return the name of the current FIO target of the current thread in the Julia environment, 
if the current thread has no target the return value will be `""`.
"""
fiotargetname() = FIO_STATE.threadtargetname[Threads.threadid()]

"""
    fiosave(object; name::String)

Save the given `object` to a file with the given `name` in the current project directory defined in `fiodir()`, 
the file name will be in the format of `{name}.dat`.

### Input
- `object` The object to be saved.
- `name` The name of the file to be saved.

### Output
The path to the file saved, in the format of `{project directory}/{name}.dat
"""
function fiosave(object; name::String)
    storageobject = object
    type = typeof(object)
    # We would like to see if there are any storage type registered for the parametric type first, 
    # for example NormalFock{Region} vs NormalFock, we will flavor the former first.
    if haskey(STORAGE_TYPES, type)
        storageobject = convert(STORAGE_TYPES[type.name.wrapper], object)
    elseif haskey(STORAGE_TYPES, type.name.wrapper)
        storageobject = convert(STORAGE_TYPES[type], object)
    end
    filepath = joinpath(fiodir(), "$name.dat")
    # Set the target name for the current thread.
    FIO_STATE.threadtargetname[Threads.threadid()] = name
    watchprogress(desc="fiosave lowering ($name)")
    jsonstring = JSON.json(storageobject)
    unwatchprogress()
    # Reset the target name for the current thread.
    FIO_STATE.threadtargetname[Threads.threadid()] = ""
    lzwcompressed::Vector = lzwcompress(jsonstring)
    open(filepath, "w") do io
        write(io, lzwcompressed)
        @info "Saved $type object to $filepath"
    end
    return filepath
end
export fiosave

"""
    fioload(name::String)

Load the object from the file with the given `name` from the current project directory defined in `fiodir()`, 
the file to be loaded will be in the path `{project directory}/{name}.dat`. Depands on the type some object might 
require extra data file to be loaded.

### Input
- `name` The name of the file to be loaded.

### Output
The object loaded from the file.
"""
function fioload(name::String)
    filepath = joinpath(fiodir(), "$name.dat")
    !isfile(filepath) && error("Object $name does not exist in $(fiodir())!")
    datasize = filesize(filepath)
    lzwcompressed = Vector{Int32}(undef, datasize/sizeof(Int32)|>round|>Integer)
    read!(filepath, lzwcompressed)
    jsonstring = lzwdecompress(lzwcompressed)
    watchprogress(desc="fioload parsing ($name)")
    object = JSON.parse(jsonstring)|>fioparse
    unwatchprogress()
    type = typeof(object)
    if haskey(PARSED_TYPES, type)
        object = convert(PARSED_TYPES[type], object)
    elseif haskey(PARSED_TYPES, type.name.wrapper)
        object = convert(PARSED_TYPES[type.name.wrapper], object)
    end
    @info "Loaded $type object from $filepath"
    return object
end
export fioload

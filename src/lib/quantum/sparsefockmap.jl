# ======================================================================================================================================
# SparseFockMap definition
"""
    SparseFockMap(outspace::FockSpace, inspace::FockSpace, rep::SparseMatrixCSC{ComplexF64, Int64})
    SparseFockMap(outspace::FockSpace, inspace::FockSpace, rep::AbstractArray{<:Number})
    SparseFockMap(outspace::FockSpace, inspace::FockSpace, mapping::Dict{Tuple{Mode, Mode}, ComplexF64})
    SparseFockMap(fockmap::FockMap; outspace::FockSpace = fockmap|>getoutspace, inspace::FockSpace = fockmap|>getinspace)

Represents an mapping between two fockspaces with the same span of the underlying Hilbert space.

### Input
- `outspace` The output `FockSpace` of this map. If this object is the multiplier, then this will be the `outspace` of the 
             resulting `FockMap`; if this object is the factor, then this must have the same span as the `inspace` of the 
             multiplier.
- `inspace`  The input `FockSpace` of this map. If this object is the multiplier, then must have the same span as the `outspace` 
             of the multiplier; if this object is the factor, then this will be the `inspace` of the resulting `FockMap`.
- `rep`      A complex sparse matrix represents the 2-point maps between the elements of the `inspace` & `outspace`.
- `mapping`  The values of the map have to be specified for a distinct 2-point pair, keyed by the pair in `Tuple`.
"""
struct SparseFockMap{A <: FockSpace, B <: FockSpace} <: FockMap{A, B}
    outspace::A
    inspace::B
    rep::SparseMatrixCSC{ComplexF64, Int64}
end
export SparseFockMap

SparseFockMap(outspace::FockSpace, inspace::FockSpace, rep::AbstractArray{<:Number}) = SparseFockMap(
    outspace, inspace, SparseMatrixCSC{ComplexF64, Int64}(rep))

function SparseFockMap(outspace::FockSpace, inspace::FockSpace, mapping::Dict{Tuple{Mode, Mode}, T})::SparseFockMap where {T <: Complex}
    rep::SparseMatrixCSC{ComplexF64, Int64} = spzeros(dimension(outspace), dimension(inspace))
    for ((out_mode::Mode, in_mode::Mode), value::ComplexF64) in mapping
        rep[outspace[out_mode], inspace[in_mode]] = value
    end
    return SparseFockMap(outspace, inspace, rep)
end

function SparseFockMap(
    fockmap::FockMap;
    outspace::FockSpace = fockmap|>getoutspace, inspace::FockSpace = fockmap|>getinspace, performpermute::Bool = true)

    mat::SparseMatrixCSC = performpermute ? permute(fockmap, outspace=outspace, inspace=inspace)|>rep : fockmap|>rep
    return SparseFockMap(outspace, inspace, mat)
end
# ======================================================================================================================================

# ==============================================================================================================================================
# FockMap default implementation
# Using SparseFockMap as the default implementation of FockMap
FockMap(outspace::FockSpace, inspace::FockSpace, rep::SparseMatrixCSC{ComplexF64, Int64}) = SparseFockMap(outspace, inspace, rep)
FockMap(outspace::FockSpace, inspace::FockSpace, rep::AbstractArray{<:Number}) = SparseFockMap(outspace, inspace, rep)
FockMap(outspace::FockSpace, inspace::FockSpace, mapping::Dict{Tuple{Mode, Mode}})::SparseFockMap= SparseFockMap(outspace, inspace, mapping)
FockMap(fockmap::FockMap; outspace::FockSpace = fockmap|>getoutspace, inspace::FockSpace = fockmap|>getinspace, performpermute::Bool = true) = (
    SparseFockMap(fockmap, outspace=outspace, inspace=inspace, performpermute=performpermute))
# ==============================================================================================================================================

# ===========================================================================================
# FockMap interface implementations
Base.:convert(::Type{SparseMatrixCSC{ComplexF64, Int64}}, source::SparseFockMap) = source.rep

getoutspace(fockmap::SparseFockMap) = fockmap.outspace

getinspace(fockmap::SparseFockMap) = fockmap.inspace
# ===========================================================================================

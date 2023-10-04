abstract type Transformation{T} <: Element{T} end
export Transformation

struct BasisFunction <: Element{Vector{Complex}}
    rep::Vector{Complex}
    dimension::Integer
    rank::Integer
end
export BasisFunction

struct AffineTransform <: Transformation{Matrix{Float64}}
    localspace::AffineSpace
    shiftvector::Vector{Float64}
    transformmatrix::Matrix{Float64}
    eigenfunctions::Dict{Tuple, BasisFunction}
    eigenvaluehashdenominator::Integer
    antiunitary::Bool
end
export AffineTransform

struct Scale <: Transformation{Matrix{Float64}}
    rep::Matrix{Float64}
end
export Scale

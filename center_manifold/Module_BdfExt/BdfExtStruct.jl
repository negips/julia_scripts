@doc raw"""
      mutable struct BdFields

      Has fields:
      Ord::Int
      N::Int
      Flds::Matrix{T}         # Real/Complex
      Coeffs::Vector{Float64}

"""
mutable struct BdFields{T<:Number}

  Ord::Int
  N::Int
  Flds::Matrix{T}
  Coeffs::Vector{Float64}

end
#---------------------------------------------------------------------- 
function BdFields(Ord::Int,v::AbstractVector{T}) where {T<:Number}

  N    = length(v)
  flds = zeros(T,N,Ord-1)
  coeffs = GetBDF(Ord)
  bdfields = BdFields{T}(Ord,N,flds,coeffs)

  return bdfields
end 
#----------------------------------------------------------------------
@doc raw"""
      mutable struct ExtFields

      Has fields:
      Ord::Int
      N::Int
      Flds::Matrix{T}         # Real/Complex
      Coeffs::Vector{Float64}

"""
mutable struct ExtFields{T<:Number}

  Ord::Int
  N::Int
  Flds::Matrix{T}
  Coeffs::Vector{Float64}

end
#---------------------------------------------------------------------- 
function ExtFields(Ord::Int,v::AbstractVector{T}) where {T<:Number}

  N         = length(v)
  flds      = zeros(T,N,Ord-1)
  coeffs    = GetEXT(Ord)
  extfields = ExtFields{T}(Ord,N,flds,coeffs)

  return extfields
end 
#----------------------------------------------------------------------






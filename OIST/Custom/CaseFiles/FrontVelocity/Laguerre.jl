# Laguerre related functions

"""
    LaguerreSpectral{T}

The spectral basis corresponding to Laguerre Polynomials with Gauss quadrature in [0,∞)
with scalar type `T`.
"""
@auto_hash_equals struct LaguerreSpectral{T} <: SpectralBasis{Line}
    nodes::Vector{T}
    weights::Vector{T}
    D::Matrix{T}

    function LaguerreSpectral(nodes::Vector{T}, weights::Vector{T}, D::Matrix{T}) where {T}
      # @argcheck length(nodes) == length(weights) == size(D,1) == size(D,2)
      new{T}(nodes, weights, D)
    end
end
#---------------------------------------------------------------------- 

"""
    LaguerreSpectral(p::Int, T=Float64)

    Generate the `Laguerre` basis of degree `p` with scalar type `T`.
"""
function LaguerreSpectral(p::Int, T=Float64)
  if p == 0
    nodes   = T[0]
    weights = T[1]        # Need to check this
  else
    nodes, weights = gausslaguerre(p+1)
  end
  D = LaguerreDerivativeMatrix(nodes,p)
  LaguerreSpectral(nodes, weights, D)
end
#---------------------------------------------------------------------- 
"""
    LaguerreDerivative(x::T, p::Int) where {T}

    Calculate the Laguerre polynomial derivative Lp'(x)
"""
function LaguerreDerivative(x::T, p::Int) where {T}

  @argcheck x >  T(0)
  @argcheck p >= 0

  if x == T(0)
    dL = -T(p) 
  else
    if p == 0
      dL = T(0)
    else
      Lp  = laguerre(x,p)
      Lp1 = laguerre(x,p-1)
      dL  = p*(Lp - Lp1)
    end  
  end

  return dL
end
#---------------------------------------------------------------------- 

"""
    LaguerreDerivativeMatrix(x::AbstractVector{T}, p::Int) where {T}

    Calculate the Derivative matrix Lp'(x) for all orders 0..p at points x
"""
function LaguerreDerivativeMatrix(x::AbstractVector{T}, p::Int) where {T}

  @argcheck minimum(x) >=  T(0)
  @argcheck p >= 0

  l = length(x)
  D = zeros(T,l,p+1)

  for j in 0:p
    for i in 1:l
      D[i,j+1] = LaguerreDerivative(x[i], j)
    end
  end

  return D
end
#---------------------------------------------------------------------- 
"""
    LaguerreTransforms(x::AbstractVector{T}, p::Int) where {T}

    Modal to Nodal (and inverse) transforms for the Lagurre Polynomials.
    Nodal values at the gauss points
"""
function LaguerreTransforms(p::Int,T=Float64)

  @argcheck p >= 0

  l   = p+1
  M2N = zeros(T,l,p+1)

  x,w = gausslaguerre(p+1)

  for j in 0:p
    for i in 1:l
      M2N[i,j+1] = laguerre(x[i], j)
    end
  end

  N2M = inv(M2N)

  return M2N, N2M
end
#---------------------------------------------------------------------- 













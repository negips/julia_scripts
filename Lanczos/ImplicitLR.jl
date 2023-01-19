function ImplicitLR!(T::Matrix,λ)
# Lower Bulge Chase Algorithm with Oblique projectors

  # H       - Tridiagonal Matrix
  # λ       - Shift

  r,c   = size(T)
  el    = eltype(T[1])
  one   = el(1)
  tol   = 1000.0*eps(abs(one))
  αtol   = 1.0e-2

  V     = Matrix{el}(I,r,c)
  W     = Matrix{el}(I,r,c)

  W,V   = SimilarityTransformBulge!(T,λ,1)
#  W[1:2,1:2] = w
#  V[1:2,1:2] = v

  Vi,Wi = ChaseBulgeTriDiagonal2!(T,λ)

# Collect Right multipliers 
  V = V*Vi
# Collect Left multipliers
  W = Wi*W

  return V,W
end
#----------------------------------------------------------------------
function ImplicitUL!(T::Matrix,λ)
  # Lower Bulge Chase Algorithm with Oblique projectors

  # H       - Tridiagonal Matrix
  # λ       - Shift

  r,c   = size(T)
  el    = eltype(T[1])
  one   = el(1)
  tol   = 1000.0*eps(abs(one))
  αtol   = 1.0e-2

  V     = Matrix{el}(I,r,c)
  W     = Matrix{el}(I,r,c)

#  W,V   = SimilarityTransformBulge!(T,λ,1)
  W,V  = CreateSuperBulgeBottomOblique(T,λ)

  Vi,Wi = ChaseBulgeUpTriDiagonal!(T,λ)

# Collect Right multipliers 
  V = V*Vi
# Collect Left multipliers
  W = Wi*W

  return V,W
end
#----------------------------------------------------------------------
function ImplicitLRSeq!(T::Matrix,μ0,nμ::Int)
# Sequential Bulge Chase

  # T       - TriDiagonal Matrix
  # nμ      - No of Shifts
  # μ       - Shifts

  r,n = size(T)

  V   = Matrix{eltype(T)}(I,n,n)
  W   = Matrix{eltype(T)}(I,n,n)

  npass = 1
# Multiple npass destroys the structure  

  if nμ == 0
    μ = [0.]
    nμ = 1
  else
    μ = μ0
  end

  for j in 1:npass
  for i in 1:nμ
    λ     = μ[i]
    Vi,Wi = ImplicitLR!(T,λ)  

#   Collect Left multipliers
    W = Wi*W  

#   Collect Right multipliers
    V = V*Vi
  end  
  end

  return V,W
end

#----------------------------------------------------------------------

function ImplicitULSeq!(T::Matrix,μ0,nμ::Int)
# Sequential Bulge Chase

  # T       - TriDiagonal Matrix
  # nμ      - No of Shifts
  # μ       - Shifts

  r,n = size(T)

  V   = Matrix{eltype(T)}(I,n,n)
  W   = Matrix{eltype(T)}(I,n,n)

  npass = 1
# Multiple npass destroys the structure  

  if nμ == 0
    μ = [0.]
    nμ = 1
  else
    μ = μ0
  end

  for j in 1:npass
  for i in 1:nμ
    λ     = μ[i]
    Vi,Wi = ImplicitUL!(T,λ)  

#   Collect Left multipliers
    W = Wi*W  

#   Collect Right multipliers
    V = V*Vi
  end  
  end

  return V,W
end

#----------------------------------------------------------------------


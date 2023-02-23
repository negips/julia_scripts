function ImplicitLR!(T::Matrix,λ,nc)
# Lower Bulge Chase Algorithm with Oblique projectors

  # H       - Tridiagonal Matrix
  # λ       - Shift
  # nc      - no of columns to go through

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

  Vi,Wi = ChaseBulgeTriDiagonalLR!(T,λ,nc)

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
    cm    = n - (i-1)
    Vi,Wi = ImplicitLR!(T,λ,cm)  

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
function ChaseBulgeTriDiagonalLR!(H::Matrix,λ,cmax)

# Modify individual entries
# I assume the the bulge size is 1
# cmax  - Maximum column number

  el = eltype(H[1])
  one = el(1)

  tol = 1000*eps(abs.(one))

  row,col = size(H)
  if row!=col
    display("H needs to be a square matrix: size(H)= $r,$c")
  end

  W  = Matrix{el}(I,row,col)
  V  = Matrix{el}(I,row,col)
  I0 = Matrix{el}(I,row,col)

  one = el(1)
  zro = el(0)

  for i in 1:cmax-2

    x1 = H[i+1,i]
    x2 = H[i+2,i]

    if abs(x1)<1.0e-4  #tol
      x1a = abs(x1)
      x2a = abs(x2)
#      println("Possible Division by 0 in Chase2 abs(x1)=$x1a, abs(x2)=$x2a, i=$i")
#      Ql,Qr = SmallX1_fix!(H,i)
#      Ql,Qr = SimilarityTransformBulge!(H,λ,i+1)     
#      V = V*Qr
#      W = Ql*W
    end

    if abs(x2)<tol
#      println("Recreating Reflector")
#      Ql,Qr = SmallX1_fix!(H,i)
#      V = V*Qr
#      W = Ql*W
    end
      
    Ql,Qr   = SimilarityTransform!(H,i,col)
    V = V*Qr
    W = Ql*W
   
  end

  return V,W
end
#----------------------------------------------------------------------



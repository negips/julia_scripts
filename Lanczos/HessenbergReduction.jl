function UpperHessenbergReduction(A::Matrix)

# Only for Square matrices 
  r,c = size(A)
  el  = eltype(A[1])
  zro = el(0)

  V   = Matrix{el}(I,r,c)
  B   = zeros(el,r,c)

  if r!=c
    display("A needs to be a square matrix: size(A)= $r,$c")
  end

  v0 = zeros(el,r)
  for i in c-2:-1:c-2 #c-2
    x             = copy(A[:,i])
    v             = HouseHolderReflector(x,i+1) 
    Q             = (I - el(2)*v*v')
    mul!(B,Q,A)
#    A            .= B
    mul!(A,B,Q')
    V[:,i]        = v 
  end

  return V
end
#----------------------------------------------------------------------
function HouseHolderReflector(x,i)

  el              = eltype(x[1])
  r               = length(x)
  v               = zeros(el,r)
  v[i]            = el(1)
  v[i:r]         .= sign(x[i])*norm(x[i:r])*v[i:r] .+ x[i:r]
  v               = v./norm(v)

  return v
end
#----------------------------------------------------------------------
function LowerHessenbergtoTriDiagonal!(H::Matrix)

  el = eltype(H[1])
  one = el(1)

  r,c = size(H)
  if r!=c
    display("H needs to be a square matrix: size(H)= $r,$c")
  end

  W  = zeros(el,r,c)
  V  = zeros(el,r,c)

  for i in 1:c-2
    w       = transpose(H[i,:])
    w[1:i]  = zeros(el,1,i)
    v       = H[:,i]
    v[1:i]  = zeros(el,i)
    v[i+1]  = one/w[i+1]
    β       = w*H[:,i]
    w       = w./β
    Q       = I - v*w
    A       = Q*H
    V[:,i+1] = -v*w[i+1]
    V[i+1,i+1] = one + V[i+1,i+1]

#   Build and store inverse Q as vectors in W
    W[i+1,i+1]  = one/V[i+1,i+1]
    for j in i+2:r
      W[j,i+1] = - V[j,i+1]/V[i+1,i+1] 
    end

    Q[:,i+1]  = W[:,i+1]
    H        .= A*Q
    
  end

  return V,W
end
#----------------------------------------------------------------------
function UpperHessenbergtoTriDiagonal!(H::Matrix)

  el = eltype(H[1])
  one = el(1)

  r,c = size(H)
  if r!=c
    display("H needs to be a square matrix: size(H)= $r,$c")
  end

  W  = zeros(el,r,c)
  V  = zeros(el,r,c)

  for i in 1:r-2
    w       = transpose(H[i,:])
    w[1:i]  = zeros(el,1,i)
    v       = H[:,i]
    v[1:i]  = zeros(el,i)
    w[i+1]  = one/v[i+1]
    β       = transpose(H[i,:])*v
    v       = v./β
    Q       = I - v*w
    A       = H*Q
    W[i+1,:] = -v[i+1]*w
    W[i+1,i+1] = one + W[i+1,i+1]

#    H      .= A 

#   Build and store inverse Q as vectors in W
    V[i+1,i+1]  = one/W[i+1,i+1]
    for j in i+2:c
      V[i+1,j] = - W[i+1,j]/W[i+1,i+1] 
    end

    Q[i+1,:]  = V[i+1,:]
    H      .= Q*A
    
  end

  return V,W
end
#----------------------------------------------------------------------
# Bulge Chase Algorithm with Oblique projectors
function NegiAlg(T::Matrix,λ)
# Also known as BulgeChase Algorithm  

  # H       - Tridiagonal Matrix
  # λ       - Shift

  r,c   = size(T)

# Create a bulge in the matrix
  T,V0,W0  = CreateLowerBulgeOblique(T,λ)
  Vi,Wi    = LowerHessenbergtoTriDiagonal!(T)

# Collect Left multipliers 
  V = Vi*V0
# Collect Right multipliers
  W = W0*Wi

  return T,V0,W0
end
#----------------------------------------------------------------------
function CreateLowerBulgeOblique(T::Matrix,λ)

  # T       - Tri-diagonal Matrix
  # λ       - Shifts

  T         = T - λ*I
  el        = eltype(T[1])
  zro       = el(0)
  one       = el(1)
  r,c       = size(T)
  x         = zeros(el,r)
  x[1]      = zro
  x[2]      = T[2,1]
  y         = transpose(T[1,:])
  α         = (y*x)[1]
  y         = y/α      # This makes <y,x> = 1.0
  β         = (y*T[:,1])[1]
  V         = I - (x*y)/β      # => <y,x/β> = 1.0/β
  Vi        = I + (x*y)/(β*(one - one/β))

# This creates the Bulge  
# A = Q0*H*Q0  
  tmp       = V*T
  T        .= tmp*Vi
  T         = T + λ*I

  return T,V,Vi
end

#----------------------------------------------------------------------
function CreateUpperBulgeOblique(T::Matrix,λ)

  # T       - Tri-diagonal Matrix
  # λ       - Shifts

  T         = T - λ*I
  el        = eltype(T[1])
  zro       = el(0)
  one       = el(1)
  r,c       = size(T)
  y         = zeros(el,1,c)
  y[1]      = zro
  y[2]      = T[1,2]
  x         = T[:,1]
  α         = (y*x)[1]
  y         = y/α      # This makes <y,x> = 1.0
  β         = (transpose(T[1,:])*x)[1]
  W         = I - (x*y)/β      # => <y,x/β> = 1.0/β
  Wi        = I + (x*y)/β/(one - one/β)

# This creates the Bulge  
# A = Q0*H*Q0  
  tmp       = T*W
  T        .= Wi*tmp
  T         = T + λ*I

  return T,Wi,W
end

#----------------------------------------------------------------------







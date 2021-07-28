# Bulge Chase Algorithm
function FrancisAlg(H::Matrix,μ::Vector,nμ::Int)
# Also known as BulgeChase Algorithm  

  # nμ      - No of Shifts
  # μ       - Shifts
  # H       - Hessenberg Matrix

  r,c   = size(H)

# Create a bulge in the matrix 
  H0,Q0 = CreateBulge(H,μ,nμ)

# Chase Bulge and return to Hessenberg form  
# Should probably code it manually  
#  hn    = hessenberg(H0)
#  Qn    = convert(Matrix,hn.Q)
#  Hn    = convert(Matrix,hn.H)      # Hn = Qn'*H0*Qn = Qn'*Q0*H*Q0'*Qn
## Collect Left multipliers 
#  Q     = Qn'*Q0                    # Hn = Q'*H*Q

  Hn,Qn  = ChaseBulgeDown(H0,nμ)
# Collect Left multipliers 
  mul!(Q,Qn,Q0)       # Q = Qn*Q0

  return Hn,Q
end
#----------------------------------------------------------------------

function FrancisSeq(H::Matrix,μ0::Vector,nμ::Int)
# Sequential Bulge Chase

  # nμ      - No of Shifts
  # μ       - Shifts
  # H       - Hessenberg Matrix

  r,n = size(H)

  Q   = Matrix{typeof(H[1,1])}(1.0I,n,n)
  Qn  = Matrix{typeof(H[1,1])}(1.0I,n,n)
  Hn  = deepcopy(H)    

  if nμ == 0
    μ = [0.]
    nμ = 1
  else
    μ = μ0
  end

  j      = 0
  nloops = 1
  βk     = 1.0
  tol    = 1.0e-12


#  println("Francis' Algorithm: nμ=$nμ, MaxLoops:$nloops, Tol=$tol")

  for i in 1:nμ
#   Create a bulge in the top left of the matrix
    λ     = μ[i:i]
    nλ    = 1
    H0,Q0 = CreateBulge(Hn,λ,nλ)

#   Collect Left multipliers 
    mul!(Qn,Q0,Q)       # Qn = Q0*Q

#   Chase Bulge through bottom right and return to 
#   Hessenberg form.  
#   Should probably code it manually  
#    hn    = hessenberg(H0)
#    Q0    = convert(Matrix,hn.Q)
#    Hn    = convert(Matrix,hn.H)      # Hn = Qn'*H0*Qn = Qn'*Q0*H*Q0'*Qn
##   Collect Left multipliers 
#    mul!(Q,Q0',Qn)       # Q = Q0'*Qn

    Hn,Q0  = ChaseBulgeDown(H0,nμ)  
##   Collect Left multipliers 
    mul!(Q,Q0,Qn)       # Q = Q0*Qn

  end
#  println("Francis Algorithm nloops: $j")

  return Hn,Q
end

#----------------------------------------------------------------------
function RevFrancisSeq(H::Matrix,μ0::Vector,nμ::Int)
# Sequential Bulge Chase

  # nμ      - No of Shifts
  # μ       - Shifts
  # H       - Hessenberg Matrix

  r,n = size(H)

  Q   = Matrix{typeof(H[1,1])}(1.0I,n,n)
  Qn  = Matrix{typeof(H[1,1])}(1.0I,n,n)
  Hn  = deepcopy(H)    

  if nμ == 0
    μ = [0.]
    nμ = 1
  else
    μ = μ0
  end

  j      = 0
  nloops = 1
  βk     = 1.0
  tol    = 1.0e-12

  println("Rev. Francis' Algorithm: nμ=$nμ, MaxLoops:$nloops, Tol=$tol")

  while βk>tol && j < nloops
    j = j+1
    for i in 1:nμ
#     Create a bulge in the top left of the matrix
      λ     = μ[i:i]
      nλ    = 1
      H0,Q0 = CreateLowerBulge(Hn,λ,nλ)

#     Collect Right multipliers 
      mul!(Qn,Q,Q0)       # Qn = Q0*Q

#     Chase Bulge through bottom right and return to 
#     Hessenberg form.  
#     Should probably code it manually  
#      hn    = hessenberg(H0)
#      Q0    = convert(Matrix,hn.Q)
#      Hn    = convert(Matrix,hn.H)      # Hn = Qn'*H0*Qn = Qn'*Q0*H*Q0'*Qn
##     Collect Left multipliers 
#      mul!(Q,Q0',Qn)       # Q = Q0'*Qn

      c1 = 2
      Hn,Q0  = ChaseBulgeUp(H0)  
#     Collect Right multipliers 
      mul!(Q,Qn,Q0)       # Q = Qn*Q0
    end
    βk = abs(Hn[r-nμ+1,r-nμ])
  end
  println("Rev. Francis Algorithm nloops: $j")

  return Hn,Q
end

#----------------------------------------------------------------------

function CreateBulge(H::Matrix,μ::Vector,nμ::Int)

  # nμ      - No of Shifts
  # μ       - Shifts
  # H       - Hessenberg Matrix

# x = p(H)*e1 = (H - μnI)...(H - μ2I)(H - μ1I)*e1
  x         = polyHe1(H,μ,nμ)
  n         = length(x)
  k1        = 1         # Zeros after this
  k2        = nμ+2      # Last non-zero entry position
  Q0,y,τ    = CreateReflectorZeros(x,k1,n)
#  Q0,y,τ    = CreateReflectorZeros2(x,k1,k2,n)

  A   = deepcopy(H)
  B   = deepcopy(H)

# This creates the Bulge  
# A = Q0*H*Q0  
  mul!(B,Q0,A)
  mul!(A,B,Q0')
 
  return A,Q0          # Return Left Multipliers
end

#----------------------------------------------------------------------
function CreateLowerBulge(H::Matrix,μ::Vector,nμ::Int)

  # nμ      - No of Shifts
  # μ       - Shifts
  # H       - Hessenberg Matrix

# x = p(H)*e1 = (H - μnI)...(H - μ2I)(H - μ1I)*e1
  r,c       = size(H)
  x         = enpolyH(H,μ,nμ)
  n         = length(x)
  Q0,y,τ    = AdjointReflectorZeros(x',r,r)

  A   = copy(H)
  B   = copy(H)

# This creates the Bulge  
# A = Q0*H*Q0  
  mul!(B,A,Q0)
  mul!(A,Q0',B)
 
  return A,Q0           # Right multiplier
end

#----------------------------------------------------------------------

#
function polyHe1(H::Matrix,μ0::Vector,nμ::Int)

# calculate x = p(H)*e1 = (H - μnI)...(H - μ2I)(H - μ1I)*e1
# when H is has a Hessenberg Structure

  # nμ      - No of Shifts
  # μ0      - Shifts
  # H       - Hessenberg Matrix

  r,c = size(H)    

  x0    = 0. *H[:,1]
  x0[1] = 1.0
  type  = typeof(H[1,1])
  x     = zeros(type,length(x0)) 

  if nμ == 0
    μ = [0.]
    nμ = 1
  else
    μ = μ0
  end  

  # calculate x = p(A)*e1 = (A - μnI)...(A-μ2I)(A-μ1I)*e1
  for i in 1:nμ
    for j in 1:i+1
      x[j] = sum(H[j,1:i].*x0[1:i]) - μ[i]*x0[j]
    end  
    x0 = 1.0*x
  end  

  return x0
end

#----------------------------------------------------------------------
function enpolyH(H::Matrix,μ0::Vector,nμ::Int)

# calculate x = en'*p(H) = en'*(H - μ1I)(H - μ2I)...(H - μnI)
# when H is has a Hessenberg Structure

  # nμ      - No of Shifts
  # μ0      - Shifts
  # H       - Hessenberg Matrix

  r,c = size(H)    

  x0    = 0. *H[:,1]'
  x0[c] = 1.0
  type  = typeof(H[1,1])
  x     = 0.0*x0 

  if nμ == 0
    μ = [0.]
    nμ = 1
  else
    μ = μ0
  end  

  # calculate x = p(A)*e1 = (A - μnI)...(A-μ2I)(A-μ1I)*e1
  for i in 1:nμ
    for j in r-i:r
      x[j] = sum(x0[j:c].*H[j:r,j]) - μ[i]*x0[j]
    end  
    x0 = 1.0*x
  end  

  return x0
end

#----------------------------------------------------------------------

function polyAe1(A::Matrix,μ0::Vector,nμ::Int)

# calculate x = p(A)*e1 = (A - μnI)...(A-μ2I)(A-μ1I)*e1
# for a general A

  # nμ      - No of Shifts
  # μ0      - Shifts
  # A       - Square Matrix

  r,c = size(A)    

  x0    = 0. *A[:,1]
  x0[1] = 1.0
  type  = typeof(A[1,1])

  if nμ == 0
    μ = [0.]
    nμ = 1
  else
    μ = μ0
  end  

  # calculate x = p(A)*e1 = (A - μnI)...(A-μ2I)(A-μ1I)*e1
  for i in 1:nμ
    x0 = A*x0 - μ[i]*x0
  end  

  return x0
end

#----------------------------------------------------------------------
function CreateReflector(x::Vector,y::Vector)

# Create Unitary matrix for reflection:
#           x --> y
# i.e.:     y = Q*x
# General Function for Real or Complex vectors 
  
  τ         = x'*y
  τ_norm    = abs(τ)
  if τ_norm>1.0e-12
    expiθ  = τ/τ_norm
  else
    println("Small τ: $τ_norm in CreateReflector")
    expiθ  = 1.0 + 0.0im
  end  

  w       = expiθ*x - y
  β       = sqrt(w'*w)
  w      .= w/β
  Q       = expiθ*(I - 2.0*w*w')

  return Q
end
#----------------------------------------------------------------------
function CreateReflectorZeros(x::Vector,k::Int,n::Int)

# Create Unitary matrix which introduces zeros after the
# kth position in the vector x
# General Function for Real or Complex vectors 
#
# x   -     Vector to Reflect
# n   -     Length of the vector
# k   -     Position after which we want zeros
#

  tol = 1.0e-12
 
  Q         = Matrix{typeof(x[1])}(1.0I,n,n)
  if k>=n
    return Q
  end

  θ         = 0.0*π/4.0

  w         = 0.0*x
  w[k]      = x[k] - norm(x[k:n])*exp(im*θ)
  w[k+1:n]  = x[k+1:n]
  β         = w'x
  τ         = 1.0/β
  Q         = (I - τ*w*w')

#  y         = Q*x
#  if norm(y[k+1:n])>tol
#    println("Remaking Projector")
#    w2         = 0.0*y
#    w2[k]      = y[k] - norm(y[k:n])*exp(im*θ)
#    w2[k+1:n]  = y[k+1:n]
#    w2         = w2
#    β          = w2'x
#    τ2         = 1.0/β
#    Q2         = (I - τ*w2*w2')
#    Q1         = copy(Q)
#    Q          = Q2*Q1
#  end  

  return Q,w,τ
end
#----------------------------------------------------------------------
function CreateReflectorZeros2(x::Vector,k1::Int,k2::Int,n::Int)

# Create Unitary matrix which introduces zeros after the
# kth position in the vector x
# General Function for Real or Complex vectors 
#
# x   -     Vector to Reflect
# n   -     Length of the vector
# k1  -     Position after which we want zeros
# k2  -     Last non-zero entry in the vector
#

  tol = 1.0e-12
 
  Q         = Matrix{typeof(x[1])}(1.0I,n,n)
  if k1>=n
    return Q
  end

  θ         = 0.0*π/4.0

  w          = zeros(typeof(x[1]),n)
  w[k1]      = x[k1] - norm(x[k1:k2])*exp(im*θ)
  w[k1+1:k2] = x[k1+1:k2]
  β          = w'x
  τ          = 1.0/β
  Q          = (I - τ*w*w')

#  y         = Q*x
#  if norm(y[k+1:n])>tol
#    println("Remaking Projector")
#    w2         = 0.0*y
#    w2[k]      = y[k] - norm(y[k:n])*exp(im*θ)
#    w2[k+1:n]  = y[k+1:n]
#    w2         = w2
#    β          = w2'x
#    τ2         = 1.0/β
#    Q2         = (I - τ*w2*w2')
#    Q1         = copy(Q)
#    Q          = Q2*Q1
#  end  

  return Q,w,τ
end
#----------------------------------------------------------------------


function TransposeReflectorZeros(x::Vector,k::Int,n::Int)

# Create Unitary matrix which introduces zeros up to the
# k-1th position in the transposed vector x
# General Function for Real or Complex vectors 
#
# x   -     Vector to Reflect
# n   -     Length of the vector
# k   -     Position after which we want zeros
 
  Q         = Matrix{typeof(x[1])}(1.0I,n,n)
  if k>=n
    return Q
  end

  θ         = 0.0*π/4.0

  w         = 0.0*x
  w[k]      = x[k] - norm(x[1:k])*exp(im*θ)
  w[1:k-1]  = x[1:k-1]
#  w         = w
  τ         = 1.0/(transpose(x)*w)
  Q         = I - τ*w*transpose(w)

  return Q,w,τ
end
#----------------------------------------------------------------------
function AdjointReflectorZeros(x::Vector,k::Int,n::Int)

# Create Unitary matrix which introduces zeros up to the
# k-1th position in the vector x
# General Function for Real or Complex vectors
# The resulting zeros are in a vector: 
#     y = x'*Q
#
# x   -     Vector to Reflect
# n   -     Length of the vector
# k   -     Position after which we want zeros
 
  Q         = Matrix{typeof(x[1])}(1.0I,n,n)
  w         = 0.0*x
  τ         = 0.0
 
  if k>n
    return Q,w,τ
  end

  θ         = 0.0*π/4.0

  w         = 0.0*x
  w[k]      = x[k] - norm(x[1:k])*exp(im*θ)
  w[1:k-1]  = x[1:k-1]
  β         = 1.0/(w'*x)
  τ         = β'    
  Q         = I - τ*w*w'

  return Q,w,τ
end
#----------------------------------------------------------------------

function ChaseBulgeDown(H0::Matrix,nμ)
#   Chase the Bulge in the Francis Algorithm

   r,c = size(H0)
   H   = deepcopy(H0)
   A   = deepcopy(H0)
   B   = deepcopy(H0)
   Q   = Matrix{typeof(H0[1,1])}(1.0I,r,c)
   T   = Matrix{typeof(H0[1,1])}(1.0I,r,c)      # tmp

   tol = 1.0e-12

   for i in 1:c-2
     x        = H[:,i]
     xnorm    = norm(x[i+2:c])

#     if xnorm>tol
       Qi,w,τ   = CreateReflectorZeros(x,i+1,r)

#       k1 = i+1
#       k2 = i+1+nμ
#       if k2>r
#         k2 = r
#       end  
#       Qi,w,τ   = CreateReflectorZeros2(x,k1,k2,r)
      
#       A = Qi*H
#       H = A*Qi'

       mul!(A,Qi,H)
       mul!(H,A,Qi')

#      y = w'*H
#      A = H .- τ*w*y
#      y = A*w
#      H = A - τ'*y*w'

#      Collect Left Multipliers      
       mul!(T,Qi,Q)
       Q = 1.0*T
#     else
#       println("$i, Improper Hessenberg: $xnorm")
#     end

   end  

   return H,Q
end
#----------------------------------------------------------------------
function ChaseBulgeUp(H0::Matrix)
#   Chase the Bulge in the Francis Algorithm

   r,c = size(H0)
   H   = deepcopy(H0)
   A   = deepcopy(H0)
   B   = deepcopy(H0)
   Q   = Matrix{typeof(H0[1,1])}(1.0I,r,c)
   T   = Matrix{typeof(H0[1,1])}(1.0I,r,c)      # tmp

   for i in c:-1:3
     x        = H[i,:]
     y        = adjoint.(x)
     Qi,w,τ   = AdjointReflectorZeros(y,i-1,c)
     
#     A = H*Qi
#     H = Qi'*A
     mul!(A,H,Qi)
     mul!(H,Qi',A)


#    Collect Right Multipliers
     mul!(T,Q,Qi) 
     Q = 1.0*T
   end  

   return H,Q
end
#----------------------------------------------------------------------








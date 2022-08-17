# Bulge Chase Algorithm
function FrancisAlg(H::Matrix,b::Int,μ::Vector,nμ::Int)
# Also known as BulgeChase Algorithm  

  # H       - Upper Hessenberg Matrix (Banded)
  # b       - Band size
  # nμ      - No of Shifts
  # μ       - Shifts

  r,c   = size(H)

# Create a bulge in the matrix
  H0,Q0 = CreateBulge(H,b,μ,nμ)

# Chase Bulge and return to Hessenberg form  
# Should probably code it manually  
#  hn    = hessenberg(H0)
#  Qn    = convert(Matrix,hn.Q)
#  Hn    = convert(Matrix,hn.H)      # Hn = Qn'*H0*Qn = Qn'*Q0*H*Q0'*Qn
## Collect Left multipliers 
#  Q     = Qn'*Q0                    # Hn = Q'*H*Q

  Hn,Qn  = ChaseBulgeDown(H0,b)
# Collect Left multipliers 
  mul!(Q,Qn,Q0)       # Q = Qn*Q0

  return Hn,Q'
end
#----------------------------------------------------------------------
function FrancisSeq(H::Matrix,b::Int,μ0::Vector,nμ::Int)
# Sequential Bulge Chase

  # H       - Hessenberg Matrix
  # b       - Band Size
  # nμ      - No of Shifts
  # μ       - Shifts

  r,n = size(H)

  Q   = Matrix{eltype(H)}(I,n,n)
  Qn  = Matrix{eltype(H)}(I,n,n)
  Hn  = deepcopy(H)    

  if nμ == 0
    μ = [0.]
    nμ = 1
  else
    μ = μ0
  end

#  println("Francis' Algorithm: nμ=$nμ, MaxLoops:$nloops, Tol=$tol")
  for i in 1:nμ
#   Create a bulge in the top left of the matrix
    λ     = μ[i:i]
    nλ    = 1
    H0,Q0 = CreateBulge(Hn,b,λ,nλ)

#   Collect Left multipliers 
    mul!(Qn,Q0,Q)       # Qn = Q0*Q
                        # H0 = Q0*Hn*Q0'

    if b==1                        
      Hn,Q0  = ChaseBulgeDown1(H0,λ,nλ)         # Hn = Q0*H0*Q0' 
#     Collect Left multipliers 
    else
      Hn,Q0  = ChaseBulgeDown(H0,b)             # Hn = Q0*H0*Q0' 
#     Collect Left multipliers 
    end  

    mul!(Q,Q0,Qn)       # Q = Q0*Qn
  end
#  println("Francis Algorithm nloops: $j")

  return Hn,Q'
end

#----------------------------------------------------------------------
function FrancisSeqExact(H::Matrix,b::Int,μ0::Vector,nμ::Int)
# Sequential Bulge Chase
# Francis Algorithm for exact Shifts
# Using Tolerance to control how many times
# Bulge Chase is done per eigenvalue shift.

  # nμ      - No of Shifts
  # μ       - Shifts
  # H       - Hessenberg Matrix

  r,n = size(H)

  Q   = Matrix{eltype(H)}(I,n,n)
  Qn  = Matrix{eltype(H)}(I,n,n)
  Hn  = deepcopy(H)    

  if nμ == 0
    μ = [0.]
    nμ = 1
  else
    μ = μ0
  end

  j      = 0
  nloops = 5 
  γ      = 1.0
  tol    = 1.0e-10

#  println("Francis' Algorithm: nμ=$nμ, MaxLoops:$nloops, Tol=$tol")
  for i in 1:nμ
#   Create a bulge in the top left of the matrix
    j = 0
    γ = 1.0
    while γ>tol && j<nloops
#    for j in 1:nloops       
      λ     = μ[i:i]
      nλ    = 1
      H0,Q0 = CreateBulge(Hn,b,λ,nλ)

#     Collect Left multipliers 
      mul!(Qn,Q0,Q)       # Qn = Q0*Q
                          # H0 = Q0*Hn*Q0'
      if b==1
        Hn,Q0  = ChaseBulgeDown1(H0,λ,nλ)         # Hn = Q0*H0*Q0' 
      else 
        Hn,Q0  = ChaseBulgeDown(H0,b)             # Hn = Q0*H0*Q0' 
      end  
#     Collect Left multipliers 
      mul!(Q,Q0,Qn)       # Q = Q0*Qn
      γ = abs(Hn[n-i+1,n-i])
      j  = j+1
    end   

  end
#  println("Francis Algorithm nloops: $j")

  return Hn,Q'
end

#----------------------------------------------------------------------
function RevFrancisSeq(H::Matrix,μ0::Vector,nμ::Int)
# Sequential Bulge Chase

  # H       - Hessenberg Matrix
  # nμ      - No of Shifts
  # μ       - Shifts

  r,n = size(H)

  Q   = Matrix{eltype(H)}(I,n,n)
  Qn  = Matrix{eltype(H)}(I,n,n)
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
  tol    = 1.0e-08

  println("Rev. Francis' Algorithm: nμ=$nμ, MaxLoops:$nloops, Tol=$tol")

  for i in 1:nμ
    j  = 0
    βk = 1.0
    while βk>tol && j < nloops
      j = j+1
#     Create a bulge in the top left of the matrix
      λ     = μ[i:i]
      nλ    = 1
      H0,Q0 = CreateLowerBulge(Hn,λ,nλ)

#     Collect Right multipliers 
      mul!(Qn,Q,Q0)       # Qn = Q*Q0

      Hn,Q0  = ChaseBulgeUp(H0)  
#     Collect Right multipliers 
      mul!(Q,Qn,Q0)       # Q = Qn*Q0
      
      βk = abs(Hn[i+1,i])
    end
  end
  println("Rev. Francis Algorithm nloops: $j")

  return Hn,Q
end

#----------------------------------------------------------------------

function CreateBulge(H::Matrix,b::Int,μ::Vector,nμ::Int)

  # H       - Hessenberg Matrix
  # b       - Block Size
  # nμ      - No of Shifts
  # μ       - Shifts

# x = p(H)*e1 = (H - μnI)...(H - μ2I)(H - μ1I)*e1
  x         = polyHe1(H,b,μ,nμ)
  n         = length(x)
  k1        = 1        # Zeros after this
  Q0,y,τ    = CreateReflectorZeros(x,k1,n)
#  Q0,y,τ    = CreateReflectorZeros2(x,k1,k2,n)

  A   = deepcopy(H)
  B   = deepcopy(H)

# This creates the Bulge  
# A = Q0*H*Q0  
  mul!(B,Q0,H)
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

function polyHe1(H::Matrix,b::Int,μ0::Vector,nμ::Int)

# calculate x = p(H)*e1 = (H - μnI)...(H - μ2I)(H - μ1I)*e1
# when H is has a Hessenberg Structure

  # H       - Hessenberg Matrix
  # b       - Block Size
  # nμ      - No of Shifts
  # μ0      - Shifts

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
    for j in 1:(i*b+1)
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

  tol = 1.0e-14
 
  Q         = Matrix{typeof(x[1])}(I,n,n)
  if k>n
    return Q
  end

  θ         = 0.0*π/4.0

  w         = 0.0*x
  w[k]      = x[k] - norm(x[k:n])*exp(im*θ)
  if k<n
    w[k+1:n]  = x[k+1:n]
  end  
  β         = w'*x
  if (abs(β)>tol)
    τ         = 1.0/β
    Q         = (I - τ*w*w')
  else
    τ         = Inf
  end

  return Q,w,τ
end
#---------------------------------------------------------------------- 
function CreateGivens(x::Vector,k::Int,n::Int)

# Create Unitary matrix which introduces zeros after the
# kth position in the vector x
# General Function for Real or Complex vectors 
#
# x   -     Vector to Reflect
# n   -     Length of the vector
# k   -     Position after which we want zeros
#

  tol = 1.0e-14
 
  Q         = Matrix{typeof(x[1])}(1.0I,n,n)
  if k>n
    return Q
  end

  θ         = 0.0*π/4.0
  if k<n
    xnorm   = abs(sqrt(x[k]'*x[k] + x[k+1]'*x[k+1]))
  
    if abs(xnorm)<tol
      cs    = 1.0
      sn    = 0.0
    else  
      cs    = x[k]/xnorm
      sn    = x[k+1]/xnorm
    end  

    Q[k,k]        = cs'
    Q[k,k+1]      = sn'
    Q[k+1,k]      = -sn
    Q[k+1,k+1]    = cs
  elseif k==n
    cs      = sign(x[k])
    sn      = 0.0 + 0.0im

    Q[k,k]  = cs
  end

  return Q,cs,sn
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
  l         = k2-k1+1
 
  Q         = Matrix{typeof(x[1])}(1.0I,n,n)
  w         = zeros(typeof(x[1]),l)
  if k1>=n
    return Q
  end

  θ         = 0.0*π/4.0

#  w          = zeros(typeof(x[1]),n)
  w[1]       = x[k1] - norm(x[k1:k2])*exp(im*θ)
  w[2:l]     = x[k1+1:k2]
  β          = w'*(x[k1:k2])
  τ          = 1.0/β
  ww         = w*w'
  Q[k1:k2,k1:k2] = I - τ*ww

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
function BandedHessenberg(H0::Matrix,b::Int)
#  Chase the Bulge in the Francis Algorithm

#  H0  - Matrix to transform into a (banded) upper Hessenberg
#  b   - Band Size.

   r,c = size(H0)
   H   = deepcopy(H0)
   A   = deepcopy(H0)
   B   = deepcopy(H0)
   Q   = Matrix{eltype(H0)}(1.0I,r,c)
   T   = Matrix{eltype(H0)}(1.0I,r,c)

   if b<1 || b>c
     println("Illegal band size b: $b")
     return H,Q
   end   

   for i in 1:c-b
     x        = H[:,i]
     Qi,w,τ   = CreateReflectorZeros(x,i+b,r)

     A = Qi*H
     H = A*(Qi')

#     Collect Left Multipliers      
     mul!(T,Qi,Q)
     Q = 1.0*T

   end  

   return H,Q
end
#----------------------------------------------------------------------

function ChaseBulgeDown(H0::Matrix,b::Int)
#   Chase the Bulge in the Francis Algorithm

   # H0     - Matrix with Bulge
   # b      - Band of the Hessenberg Matrix to convert back to

   r,c = size(H0)
   H   = deepcopy(H0)
   A   = deepcopy(H0)
   B   = deepcopy(H0)
   Q   = Matrix{eltype(H0)}(I,r,c)
   T   = Matrix{eltype(H0)}(I,r,c)      # tmp

   tol = 1.0e-12

   for i in 1:c-b
     x       = H[:,i]

     k1      = i+b
     Qi,w,τ  = CreateReflectorZeros(x,k1,r)
     
#     k1    = i+1
#     Qi    = CreateGivens(x,k1,c)

     A = Qi*H
     H = A*(Qi')

#    Collect Left Multipliers      
     mul!(T,Qi,Q)
     Q = 1.0*T

   end  

   return H,Q
end
#----------------------------------------------------------------------
function ChaseBulgeDown1(H0::Matrix,μ,nμ::Int)
#   Chase the Bulge in the Francis Algorithm
#   Here I am assuming we have a bulge size of 1
#   i.e. nμ=1 was used for creating the Bulge and band size = 1

   r,c = size(H0)
   H   = deepcopy(H0)
   A   = deepcopy(H0)
   B   = deepcopy(H0)
   Q   = Matrix{eltype(H0)}(I,r,c)
   T   = Matrix{eltype(H0)}(I,r,c)      # tmp

   tol = 1.0e-12

   for i in 1:c-1
     x        = H[:,i]

     if i<=c-2
       xb       = x[i+2]
     else
       xb       = 1.0
     end

     if abs(xb)>tol

       k1       = i+1
       Qi,cs,sn = CreateGivens(x,k1,r)
       if k1<r 
         v1 = H[k1,:]
         v2 = H[k1+1,:]
         H[k1,:]   =  cs'*v1 + sn'*v2
         H[k1+1,:] = -sn*v1  + cs*v2

         v3 = H[:,k1]
         v4 = H[:,k1+1]
         H[:,k1]   =  v3*cs  + v4*sn
         H[:,k1+1] = -v3*sn' + v4*cs'


         v1 = Q[k1,:]
         v2 = Q[k1+1,:]
         Q[k1,:]   =  cs'*v1 + sn'*v2
         Q[k1+1,:] = -sn*v1  + cs*v2
       else
         v1 = H[k1,:]
         H[k1,:]   =  cs'*v1

         v3 = H[:,k1]
         H[:,k1]   =  v3*cs


         v1 = Q[k1,:]
         Q[k1,:]   =  cs'*v1
       end

     else  

       if i<c-2
#        Create Bulge in the next column

#        println("Invariant subspace reached: ($(i+1),$i): $xb")
#        println("Generating new Lower Bulge for i=$(i+1)")
         j = i+1
         hs,qs = CreateBulge(H[j:r,j:c],1,μ,nμ)

         Qi = Matrix{eltype(H)}(I,r,c)
         Qi[j:r,j:c] = qs

         mul!(A,Qi,H)
         mul!(H,A,Qi')
#         Collect Left Multipliers      
         mul!(T,Qi,Q)
         Q = copy(T)

       end
     end  

   end      # for i in 1:c-1 

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
     
     mul!(A,H,Qi)
     mul!(H,Qi',A)

#    Collect Right Multipliers
     mul!(T,Q,Qi) 
     Q = 1.0*T
   end  

   return H,Q
end
#----------------------------------------------------------------------








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

  Hn,Qn  = ChaseBulgeDown(H0)
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

  println("Francis' Algorithm: nμ=$nμ")

  j      = 0
  nloops = 1
  βk     = 1.0
  tol    = 1.0e-10

  while βk>tol && j < nloops
    j = j+1
    for i in 1:nμ
#     Create a bulge in the top left of the matrix
      λ     = μ[i:i]
      nλ    = 1
      H0,Q0 = CreateBulge(Hn,λ,nλ)

#     Collect Left multipliers 
      mul!(Qn,Q0,Q)       # Qn = Q0*Q

#     Chase Bulge through bottom right and return to 
#     Hessenberg form.  
#     Should probably code it manually  
#      hn    = hessenberg(H0)
#      Q0    = convert(Matrix,hn.Q)
#      Hn    = convert(Matrix,hn.H)      # Hn = Qn'*H0*Qn = Qn'*Q0*H*Q0'*Qn
##     Collect Left multipliers 
#      mul!(Q,Q0',Qn)       # Q = Q0'*Qn

      Hn,Q0  = ChaseBulgeDown(H0)  
#     Collect Left multipliers 
      mul!(Q,Q0,Qn)       # Q = Q0*Qn
    end
    βk = abs(Hn[r-nμ+1,r-nμ])
  end
  println("Francis Algorithm nloops: $j")

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

  println("Francis' Algorithm: nμ=$nμ")

  j      = 0
  nloops = 1
  βk     = 1.0
  tol    = 1.0e-10

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

      Hn,Q0  = ChaseBulgeUp(H0)  
#     Collect Right multipliers 
      mul!(Q,Qn,Q0)       # Q = Qn*Q0
    end
    βk = abs(Hn[r-nμ+1,r-nμ])
  end
  println("Francis Algorithm nloops: $j")

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
  Q0,y,τ    = CreateReflectorZeros(x,1,n)

  A   = copy(H)
  B   = copy(H)

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
    x0 = x
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
    x0 = x
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
  w         = w
  β         = w'x
  if (abs(β)>tol)
    τ       = 1.0/β
    Q       = (I - τ*w*w')
  else
    w       = 0.0*x
    τ       = 0.0
  end  

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

function ChaseBulgeDown(H0::Matrix)
#   Chase the Bulge in the Francis Algorithm

   r,c = size(H0)
   H   = copy(H0)
   A   = copy(H0)
   B   = copy(H0)
   Q   = Matrix{typeof(H0[1,1])}(1.0I,r,c)
   T   = Matrix{typeof(H0[1,1])}(1.0I,r,c)      # tmp

   for i in 1:c-2
     x        = H[:,i]
     Qi,y,τ   = CreateReflectorZeros(x,i+1,r)
     
     A = Qi*H
     H = A*Qi'

#    Collect Left Multipliers      
     T = Qi*Q
     Q = copy(T)
   end  

   return H,Q
end
#----------------------------------------------------------------------
function ChaseBulgeUp(H0::Matrix)
#   Chase the Bulge in the Francis Algorithm

   r,c = size(H0)
   H   = copy(H0)
   A   = copy(H0)
   B   = copy(H0)
   Q   = Matrix{typeof(H0[1,1])}(1.0I,r,c)
   T   = Matrix{typeof(H0[1,1])}(1.0I,r,c)      # tmp

   for i in c:-1:3
     x        = H[i,:]
     y        = adjoint.(x)
     Qi,w,τ   = AdjointReflectorZeros(y,i-1,c)
     
     A = H*Qi
     H = Qi'*A

#    Collect Right Multipliers      
     T = Q*Qi
     Q = copy(T)
   end  

   return H,Q
end
#----------------------------------------------------------------------








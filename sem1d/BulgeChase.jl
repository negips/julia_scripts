# Bulge Chase Algorithm
function FrancisAlg(H::Matrix,μ::Vector,nμ::Int)
# Also known as BulgeChase Algorithm  

  # nμ      - No of Shifts
  # μ       - Shifts
  # H       - Hessenberg Matrix

# Create a bulge in the matrix 
  H0,Q0 = CreateBulge(H,μ,nμ)

# Chase Bulge and return to Hessenberg form  
# Should probably code it manually  
  hn    = hessenberg(H0)
  Qn    = convert(Matrix,hn.Q)
  Hn    = convert(Matrix,hn.H)      # Hn = Qn'*H0*Qn = Qn'*Q0'*H*Q0*Qn

# Collect Right multipliers 
  Q     = Q0*Qn                     # Hn = Q'*H*Q
 
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

  for i in 1:nμ
#   Create a bulge in the matrix
    λ     = μ[i:i]
    nλ    = 1
    H0,Q0 = CreateBulge(Hn,λ,nλ)

#   Collect Right multipliers 
    mul!(Qn,Q,Q0)       # Qn = Q*Q0

#   Chase Bulge and return to Hessenberg form  
#   Should probably code it manually  
    hn    = hessenberg(H0)
    Q0    = convert(Matrix,hn.Q)
    Hn    = convert(Matrix,hn.H)      # Hn = Qn'*H0*Qn = Qn'*Q0'*H*Q0*Qn

#    Hn,Q0  = ChaseBulge(H0)  

#   Collect Right multipliers 
    mul!(Q,Qn,Q0)       # Q = Qn*Q0
  end

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

  A   = deepcopy(H)
  B   = deepcopy(H)

# This creates the Bulge  
# A = Q0*H*Q0  
  mul!(B,Q0,A)
  mul!(A,B,Q0')
 
  return A,Q0
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
 
  Q         = Matrix{typeof(x[1])}(1.0I,n,n)
  if k>=n
    return Q
  end

  θ         = 6.0*π/4.0

  w         = 0.0*x
  w[k]      = x[k] - norm(x[k:n])*exp(im*θ)
  w[k+1:n]  = x[k+1:n]
  w         = w
  τ         = 1.0/(w'*x)
  Q         = (I - τ*w*w')

  return Q,w,τ
end
#----------------------------------------------------------------------

function ChaseBulge(H0)
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
     v = y'*H
     A = y*v
     
     H .= H .- τ*A

#     mul!(A,Qi,H) 
#     mul!(H,A,Qi)

#    Collect Right Multipliers      
     mul!(T,Q,Qi')
     Q = copy(T)
   end  

   return H,Q
end
#----------------------------------------------------------------------








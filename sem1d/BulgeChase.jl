# Bulge Chase Algorithm
function BulgeChase(H::Matrix,μ::Vector,nμ::Int)

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

function BulgeChaseSeq(H::Matrix,μ0::Vector,nμ::Int)
# Sequential Bulge Chase

  # nμ      - No of Shifts
  # μ       - Shifts
  # H       - Hessenberg Matrix

  Q   = Matrix{typeof(A[1,1])}(1.0I,n,n)
  Qn  = Matrix{typeof(A[1,1])}(1.0I,n,n)
  Hn  = copy(H)    

  if nμ == 0
    μ = [0.]
    nμ = 1
  else
    μ = μ0
  end  

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
  x     = polyHe1(H,μ,nμ)
  xnorm = sqrt(x'*x)
  y     = 0.0*x

  y[1]      = xnorm
  Q0        = CreateReflector(x,y)    

  A   = copy(H)
  B   = copy(H)

# This creates the Bulge  
# A = Q0*H*Q0  
  mul!(B,Q0',A)
  mul!(A,B,Q0)
 
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
    expiθ  = 1.0 + 0.0im
  end  

  w       = expiθ*x - y
  β       = sqrt(w'*w)
  w      .= w/β
  Q       = expiθ*(I - 2.0*w*w')

  return Q
end




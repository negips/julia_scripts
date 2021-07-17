# Testing the Francis Algorithm

include("BulgeChase.jl")
include("ExplicitShiftedQR.jl")
using LinearAlgebra
using Random
using PyPlot

rng = MersenneTwister(1234)

close("all")

n = 8

λs = randn(rng,ComplexF64,n)

AR = randn(rng,ComplexF64,n,n)
A  = inv(AR)*diagm(λs)*AR
h = hessenberg(A)
H = convert(Matrix,h.H)

F = eigen(H)
ind = sortperm(abs.(F.values))
λ   = F.values[ind]

m   = 1

μ  = λ[1:m]
#μ  = [0.]
nμ = length(μ)

Q  = Matrix{typeof(A[1,1])}(1.0I,n,n)

Hn = copy(H)

niter = 1
res   = zeros(Float64,niter)     # Vector
res1  = zeros(Float64,niter)     # Vector

#for i in 1:niter
#  global Q, Hn, res

  i = 1

  Hpert,Q0 = CreateBulge(Hn,μ,nμ)    # = Q0'*p(H)*Q0

# Collect Right multipliers 
  Q     = Q*Q0
  
  hn    = hessenberg(Hpert)
  Qn    = convert(Matrix,hn.Q)
  Hn    = convert(Matrix,hn.H)      # Hn = Qn'*Hpert*Qn = Qn'*Q0'*H*Q0*Qn

# Collect Right multipliers 
  Q     = Q*Qn                      # Hn = Q'*H*Q


  res[i]  = abs(Hn[n-m+1,n-m])

# Test QR

  B   = H - μ[1]*I
  t   = qr(B)
  B1  = t.R*t.Q + μ[1]*I


# Hessenberg Decomposition is uniquely determined by the first column of Q
# So I am finding the reflector for the first column so that it gives the same Q as the QR decomposition
# of p(A) = QR
  QR_Q1 = t.Q[:,1]

  u = Q[:,1] - QR_Q1
  τ = u'*Q[:,1]
  Q2 = I - (1.0/τ)*u*u'

  yy = Q2*Q[:,1]
  Q3 = Q2*Q

  R       = (Q3)'*(H-μ[1]*I)
  Hnn     = R*(Q3) + μ[1]*I


  ν       = λ[1:3]
  nν      = length(ν)
  Q4,H4   = ShiftedQR(H,ν,nν,2)    

  H4[:,5:8]

#  ν       = λ[2]
#  nν      = 1
#  Q5,H5,R5 = ShiftedQR(H4,ν,nν,1)    


#  t.Q

#end

#
#semilogy(res)
#semilogy(res1)

# Hn[n-m:n,n-m:n]

#x = polyHe1(H,μ,nμ)
# y = 0.0*x
# y[1] = -sqrt(x'*x)
# 
# u = x-y
# #u = u/u[1]
# τ = u'*x
# 
# Q = I - (1.0/τ)*u*u'

#Q0*x

# x = 0.0*x
# x[1] = im
# x[2] = 1.0
# xr   = real.(x)
# xi   = imag.(x)
# 
# 
# u = y-x
# u = u/norm(u)
# Q = I - 2.0*u*u'
# 
# Q*x



# Testing the Francis Algorithm

include("BulgeChase.jl")
include("ExplicitShiftedQR.jl")
include("FrancisAlg.jl")

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

#  Hpert,Q0 = CreateBulge(Hn,μ,nμ)    # = Q0'*p(H)*Q0
  μ         = λ[1:1]
  nμ        = 1
  H1,Q1     = BulgeChase(H,μ,nμ)
  μ         = λ[2:2]
  nμ        = 1
  H2,Q2     = BulgeChase(H1,μ,nμ)
  μ         = λ[3:3]
  nμ        = 1
  H3,Q3     = BulgeChase(H2,μ,nμ)

  H22,Q22   = BulgeChaseSeq(H,λ[1:3],3)

## Collect Right multipliers 
#  Q     = Q*Q0
#  
#  hn    = hessenberg(Hpert)
#  Qn    = convert(Matrix,hn.Q)
#  Hn    = convert(Matrix,hn.H)      # Hn = Qn'*Hpert*Qn = Qn'*Q0'*H*Q0*Qn
#
## Collect Right multipliers 
#  Q     = Q*Qn                      # Hn = Q'*H*Q


  res[i]  = abs(Hn[n-m+1,n-m])

# Test QR

  B   = H - μ[1]*I
  t   = qr(B)
  B1  = t.R*t.Q + μ[1]*I


x = polyHe1(H,μ,nμ)
y = 0.0*x
xnorm = sqrt(x'*x)
y[1] = abs(xnorm)
τ    = x'*y
τ_norm = abs(τ)
if τ_norm>1.0e-12
  expiθ  = τ/τ_norm
else 
  expiθ  = 1.0 + 0.0im
end  

u    = (x-y)
unorm = sqrt(u'*u)
u   .= u/unorm
Qi   = I - 2.0*u*u'

w    = expiθ*x - y
unorm = sqrt(w'*w)
w   .= w/unorm
U    = expiθ*(I - 2.0*w*w')

U*x

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



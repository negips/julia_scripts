# Testing Eigenvalue algorithm implementation
println("Testing BiOrthogonal QR method")

using LinearAlgebra
using Random
#using Pseudospectra
using Printf
using PyPlot

include("ExplicitShiftedQR.jl")
include("BulgeChase.jl")
include("HessenbergReduction.jl")
include("NegiEig.jl")

close("all")

vt = ComplexF64
#vt = Float64

zro   = vt(0.0)

rng = MersenneTwister(1254)

n     = 200
A     = randn(rng,vt,n,n)

AH    = copy(A)
q     = HessenbergReduction!(AH)

#AT    = copy(AH)
#v0,w0 = UpperHessenbergtoTriDiagonal2!(AT)

d0    = randn(rng,vt,n)
du    = randn(rng,vt,n-1)
dl    = randn(rng,vt,n-1)
AT    = Matrix(Tridiagonal(dl,d0,du))
A     = copy(AT)


B     = copy(AT)
C     = copy(AT)
θ     = eigvals(AT)


niter = 5 
ern   = zeros(Float64,niter+1)
erqr  = zeros(Float64,niter+1)

ern[1]  = abs(AT[n,n-1])
erqr[1] = abs(C[n,n-1])

for i in 1:niter
  local λ,μ
  global v,w
  global AT,C
  global ern,erqr

  λ         = AT[n,n]
  AT,v,w    = NegiAlg2(AT,λ)
  er        = AT[n,n-1]
  l         = AT[n,n]
  ern[i+1]  = abs(er)

  μ         = zeros(vt,1)
  μ[1]      = C[n,n]
  C,Q       = FrancisAlg(C,1,μ,1) 

  erb       = C[n,n-1]
  lb        = C[n,n]
  erqr[i+1] = abs(erb) 
 
#  println("$er,   $l, $erb,   $lb")
end  

display([eigvals(B) eigvals(AT) eigvals(C)])

semilogy(ern)
semilogy(erqr)

λ = B[1,1]
b = copy(B) - λ*I
#c,wi,w = CreateUpperBulgeOblique(b,λ)
#d,x,y  = CreateUpperBulgeRightOblique(b,λ)

c,x,y   = CreateLowerBulgeOblique(b,λ)
α       = norm(x[2]*y)
y       = y/α
#x       = x./(10*abs(α))
q       = I - x*y
qi      = inv(q)
d       = q*b*qi

d       = copy(c)
e       = copy(c)

#v1,w1     = ChaseBulgeTriDiagonal!(d)
#
#v2,w2     = ChaseBulgeTriDiagonal2!(e)

#V,W     = LowerHessenbergtoTriDiagonal2!(c)
#norm(W*V - I)
#spy(d,precision=1.0e-12)



println("Done.")














# Testing Eigenvalue algorithm implementation
println("Testing Triangular QR method")

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

n     = 100
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


niter = 50
ern   = zeros(Float64,niter+1)
erqr  = zeros(Float64,niter+1)

ern[1]  = abs(AT[n,n-1])
erqr[1] = abs(C[n,n-1])

nconv   = 0

for i in 1:niter
  local λ,μ
  global v,w
  global AT,C
  global ern,erqr
  global nconv

  j = n - nconv
#  if (mod(i,2)==0)
    λ         = AT[j,j]
    AT,v,w    = NegiAlg2(AT,λ)
    er        = AT[j,j-1]
#  else
#    λ         = AT[1,1]
#    AT,v,w    = NegiAlg4(AT,λ)
#    er        = AT[2,1]
#  end
  l         = AT[j,j]
  ern[i+1]  = abs(er)

  μ         = zeros(vt,1)
  μ[1]      = C[j,j]
  C,Q       = FrancisAlg(C,1,μ,1) 

  erb       = C[j,j-1]
  lb        = C[j,j]
  erqr[i+1] = abs(erb)

  if (abs(ern[i+1])<1.0e-12)
    nconv = nconv + 1
  end  
 
#  println("$er,   $l, $erb,   $lb")
end  

display([eigvals(B) eigvals(AT) eigvals(C)])

semilogy(ern)
semilogy(erqr)

λ = zro #B[n,n]
b = copy(B) - λ*I
#c,wi,w = CreateUpperBulgeOblique(b,λ)
#d,x,y  = CreateUpperBulgeRightOblique(b,λ)

c,x1,y1   = CreateLowerBulgeOblique(b,λ)
v1,w1   = ChaseBulgeTriDiagonal2!(c)

d,x2,y2 = CreateUpperBulgeOblique(b,λ)
v2,w2   = ChaseBulgeTriDiagonal3!(d)

e        = copy(b)
e[n,n-2] = rand(vt)
f        = copy(e)

g        = copy(b)
t,x,y    = CreateLowerRightBulgeOblique(g,λ)
u        = copy(t)
v3,w3    = ChaseBulgeTriDiagonal4!(u)



#v1,w1     = ChaseBulgeTriDiagonal!(d)
#
#v2,w2     = ChaseBulgeTriDiagonal2!(e)

#V,W     = LowerHessenbergtoTriDiagonal2!(c)
#norm(W*V - I)
#spy(d,precision=1.0e-12)



println("Done.")














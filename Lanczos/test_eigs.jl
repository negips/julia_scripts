# Testing Eigenvalue algorithm implementation
println("Testing TriDiagonal QR method")

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

n     = 500
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


niter = 20
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

  j = n #- nconv
  λ         = AT[j,j]
  AT,v,w    = NegiAlg2(AT,λ)
  er        = AT[j,j-1]
  l         = AT[j,j]
  ern[i+1]  = abs(er)

  μ         = zeros(vt,1)
  μ[1]      = C[j,j]
  C,Q       = FrancisAlg(C,1,μ,1) 

  erb       = C[j,j-1]
  lb        = C[j,j]
  erqr[i+1] = abs(erb)

  if (abs(ern[i+1])<1.0e-12)
#    nconv = nconv + 1
#    break
  end

  superd = diag(AT,1)
  subd   = diag(AT,-1)

  plot(abs.(superd), linestyle="-")
  plot(abs.(subd), linestyle="--")

  pause(0.000001)
#  println("$er,   $l, $erb,   $lb")
end  

#display([eigvals(B) eigvals(AT) eigvals(C)])

figure(num=2)
semilogy(ern)
#semilogy(erqr)

#λ = zro #B[n,n]
#b = copy(B) - λ*I
#
#g        = copy(b)
##g[n,n-1] = zro
##T,x,y    = CreateUpperRightBulgeOblique(g,zro)
#T,x,y    = CreateLowerBulgeOblique(g,zro)
#
#G        = copy(T)
#ql,qr    = SimilarityTransform!(G,1,n)
#
#v1,w1     = ChaseBulgeTriDiagonal!(d)




println("Done.")














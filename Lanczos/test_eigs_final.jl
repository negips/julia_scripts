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
include("ImplicitLR.jl")
include("BiOrthoIRst.jl")

close("all")

vt = ComplexF64
#vt = Float64

zro   = vt(0.0)

rng = MersenneTwister(1254)

n     = 20

d0    = randn(rng,vt,n)
du    = randn(rng,vt,n-1)
dl    = randn(rng,vt,n-1)
AT    = Matrix(Tridiagonal(dl,d0,du))
A     = copy(AT)


B     = copy(AT)
C     = copy(AT)
θ     = eigvals(AT)


niter = 2
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
#  AT,v,w    = NegiAlg2(AT,λ)
  v,w       = ImplicitLR!(AT,λ)

  er        = AT[j,j-1]
  l         = AT[j,j]
  ern[i+1]  = abs(er)

  if (abs(ern[i+1])<1.0e-12)
#    nconv = nconv + 1
#    break
  end

#  pause(0.000001)
#  println("$er,   $l, $erb,   $lb")
end  

#display([eigvals(B) eigvals(AT) eigvals(C)])

figure(num=2)
semilogy(ern)
#semilogy(erqr)

μ,nμ = GetLowerShifts(AT,3);

vi,wi = ImplicitLRSeq!(AT,μ,nμ)
#vi,wi = ImplicitLRSeq!(AT,μ,nμ)


println("Done.")














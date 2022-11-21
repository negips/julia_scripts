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
one   = vt(1.0)

tol   = 10*eps(abs(one))

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

niter = 400
ern   = zeros(Float64,niter+1) .+ 1.0e-26
erqr  = zeros(Float64,niter+1)
mnorm = zeros(Float64,niter+1)

ern[1]  = abs(AT[n,n-1])
erqr[1] = abs(B[n,n-1])
mnorm[1]= norm(AT)

nconv   = 0

for i in 1:niter
  local λ,μ
  global v,w
  global AT,C
  global ern,erqr,mnorm
  global nconv

  j = n - nconv
  λ         = AT[j,j]
  AT[1:j,1:j],v,w    = NegiAlg2(AT[1:j,1:j],λ)
  er        = AT[j,j-1]
  l         = AT[j,j]
  ern[i+1]  = abs(er)

  mnorm[i+1] = norm(AT)
#  j = n - nconv
#  λ         = AT[j,j]
#  AT,v,w    = NegiAlg2(AT,λ)
#  er        = AT[j,j-1]
#  ern[i+1]  = abs(er)

  if abs(er)<tol
    nconv = nconv + 1

    if nconv == n-2
      break
    end  
  end  

#  println("$er,   $l, $erb,   $lb")
end  

#display([eigvals(B) eigvals(AT) eigvals(C)])

figure(num=1)
semilogy(ern)
#semilogy(erqr)

d = diag(AT,-1)
figure(num=2)
semilogy(abs.(d))

figure(num=3)
plot(mnorm,linestyle="none",marker="o")


δ = eigvals(B) .- eigvals(AT)

display(δ)
normδ = norm(δ)
println("Eigenvalues error norm = $normδ")
println("Done.")














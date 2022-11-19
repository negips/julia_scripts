# Testing IRAM implementation
println("Testing BiOrthogonal QR method")

using LinearAlgebra
using Random
#using Pseudospectra
using Printf
using PyPlot

include("ExplicitShiftedQR.jl")
include("BulgeChase.jl")
include("HessenbergReduction.jl")

close("all")

rng = MersenneTwister(1234)

vt = ComplexF64
#vt = Float64

zro   = vt(0.0)

rng = MersenneTwister(1254)

n     = 10
A     = randn(rng,vt,n,n)
q     = HessenbergReduction!(A)
v0,w0 = UpperHessenbergtoTriDiagonal!(A)

B     = copy(A)
C     = copy(A)
θ     = eigvals(A)

for i in 1:3
  local λ,μ
  global v,w
  global A,C
  λ     = A[n,n]
  A,v,w = NegiAlg(A,λ)
  er    = abs(A[n,n-1])
  l     = A[n,n]

  μ     = zeros(vt,1)
  μ[1]  = C[n,n]
  C,Q   = FrancisAlg(C,1,μ,1) 

  erb   = abs(C[n,n-1])
  lb    = C[n,n]
 
  println("$er,   $l, $erb,   $lb")
end  

#C = copy(A[1:n-1,1:n-1])
#
#m = n-1
#for i in 1:4
#  local λ
#  local v,w
#  global C
#  λ     = C[m,m]
#  C,v,w = NegiAlg(C,λ)
#  er    = C[m,m-1]
#  l     = C[m,m]
#  println("$er,   $l")
#end  

display([eigvals(B) eigvals(A)])

println("Done.")














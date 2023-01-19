# Testing IRAM implementation
println("Testing Implicit UR method")

using LinearAlgebra
using Random
#using Pseudospectra
using Printf
using PyPlot

include("BiOrthoUpd.jl")
include("BiOrthoIRst.jl")
include("ExplicitShiftedQR.jl")
include("IRBiOrtho.jl")
include("BulgeChase.jl")
include("HessenbergReduction.jl")
include("NegiEig.jl")
include("ImplicitLR.jl")
include("SetZero.jl")

close("all")

rng = MersenneTwister(1237)

vt = ComplexF64
#vt = Float64

zro   = vt(0.0)
one   = vt(1.0)

n = 8     # Matrix size

dl  = randn(vt,n-1)
du  = randn(vt,n-1)
d   = randn(vt,n)
M   = Matrix(Tridiagonal(dl,d,du))
Mnorm = norm(M)

μ   = eigvals(M);

A   = copy(M)

nλ  = 5
λ   = μ[1:nλ]

nθ  = n - nλ
θ   = μ[nλ+1:n]

v1,w1 = ImplicitLRSeq!(A,λ,nλ)
Anorm = norm(A)

B     = copy(A)
v2,w2 = ImplicitULSeq!(B,θ,nθ)
Bnorm = norm(B)

Q     = zeros(vt,n,n)
for i in 1:n
  Q[i,n-i+1] = one
end  

C     = Q*A*Q
v3,w3 = ImplicitLRSeq!(C,λ,nλ)
Cnorm = norm(C)

println("Norm(M)=$Mnorm; Norm(A)=$Anorm; Norm(B)=$Bnorm; Norm(C)=$Cnorm")
println("Done")








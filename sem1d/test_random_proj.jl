# Testing the Francis Algorithm

include("BulgeChase.jl")
include("ExplicitShiftedQR.jl")
include("BulgeChase.jl")

using LinearAlgebra
using Random
using PyPlot

close("all")

include("sem_init_ref.jl")
include("sem_main.jl")


r,c = size(OPg)

# No of random vectors
n   = 100

V   = randn(ComplexF64,r,n)
q_r = qr(V)
Q   = convert(Matrix,q_r.Q)

Hr  = Q'*OPg*Q

F   = eigen(Hr)
λ   = im*F.values

h   = figure(num=1,figsize=[8.0,6.0])
ax1 = gca()
pλ = ax1.plot(real.(λ),imag.(λ), linestyle="none",marker=".", markersize=8)


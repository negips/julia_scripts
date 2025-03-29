#!/bin/julia

println("Testing SVD")

using LinearAlgebra
using PyPlot

include("yinyun_matrices.jl")

nr    = 3
nc    = 4
n     = nr*nc
lafs  = 16

α     = 1.53e8
β     = 580.0
kai   = 2.94e7
δ     = 130.0
γ     = 1.0e2
Ca    = 10.0e-6


L     = YinyunMatrix(nr,nc,α,β,kai,δ,γ,Ca)
v0    = zeros(Float64,n)
v0[1] = 1.0
v0    = v0/norm(v0)
un    = zeros(Float64,n)
un[n] = 1.0
un    = un/norm(un)
P     = v0*v0'
Q     = un*un'

nT    = 5000
T     = LinRange(0.00001,0.5,nT)
svals = zeros(Float64,nT)

for i in 1:nT
  global svals

  t     = T[i]
  ExpLT = exp(L*t)
  A     = L*ExpLT
  QAP   = Q*A*P

  F     = svd(QAP)
  svals[i] = F.S[1]
end

close("all")
plot(T,svals)
xlabel("Time",fontsize=lafs)
ylabel(L"σ_{1}",fontsize=lafs)

println("Done")









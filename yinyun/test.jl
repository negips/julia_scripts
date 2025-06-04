#!/bin/julia

println("Testing SVD")

using LinearAlgebra
include("yinyun_matrices.jl")

nr    = 3
nc    = 4
n     = nr*nc

α     = 1.53e8
β     = 580.0
kai   = 2.94e7
δ     = 130.0
γ     = 1.0e2
Ca    = 10.0e-6

T     = 0.001
L     = YinyunMatrix(nr,nc,α,β,kai,δ,γ,Ca)
ExpLT = exp(L*T)
A     = L*ExpLT

v0    = zeros(Float64,n)
v0[1] = 1.0
v0    = v0/norm(v0)

P     = v0*v0'

un    = zeros(Float64,n)
un[n] = 1.0
un    = un/norm(un)

Q     = un*un'

QAP   = Q*A*P


F     = svd(QAP)











# Testing the module

include("SEM1D.jl")
using .SEM1D

using LinearAlgebra
using SparseArrays
using Printf
using PolynomialBases

include("AssembleMatrixGinzburgLandau.jl")


N     = 8
Nd    = 12
nel   = 4
xs    = 0.0
xe    = 15.0

# Input parameters
Inp   = SEM_Input(N,Nd,nel,xs,xe)
# Nodal Bases
B0    = LobattoLegendre(Inp.N)
Bd    = LobattoLegendre(Inp.Nd)
# Geometric Matrices
GeoM  = SEM_GeoMat(B0,Bd,Inp)
δ     = ones(ComplexF64,4)
A, B, OP, Conv, Src, Lap = AssembleMatrixGLSparse(δ,Inp,GeoM,B0)
println("Done.")

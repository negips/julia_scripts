# Testing the module

include("../Module_SEM1D/SEM1D.jl")

using LinearAlgebra
using Printf
using PolynomialBases

include("GL_Functions.jl")
#-------------------------------------------------- 

# Input parameters
Inp   = Get_SEM1D_Input()
# Nodal Bases
B0    = LobattoLegendre(Inp.N)
Bd    = LobattoLegendre(Inp.Nd)
# Geometric Matrices
GeoM  = SEM1D.SEMGeoMat(B0,Bd,Inp)

println("SEM Setup Done.")
















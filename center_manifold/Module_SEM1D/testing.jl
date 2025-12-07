# Testing the module

include("SEM1D.jl")
using .SEM1D

using LinearAlgebra
using SparseArrays
using Printf
using PolynomialBases

#include("AssembleMatrixGinzburgLandau.jl")


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

δ     = ones(ComplexF64,4)    # Parameters
ifconj = false
#if (ifconj)
#  U    = conj(U)        # -δ_1
#  γ    = conj(γ)        # δ_4
#  μ0   = conj(μ0)       # δ_2
#  μx   = conj(μx)       # δ_3
#end  

# GinzburgLandau Linear Operators
L,  B, OP,  Conv,   Src, Lap  = GinzburgLandauSparse(δ,Inp,GeoM,B0)
AL, B, AOP, AConv, ASrc, ALap = AdjointGinzburgLandauSparse(δ,Inp,GeoM,B0)

ifperiodic = false
ndof, glnum = SEM_Global_Num(GeoM.xm1,ifperiodic)
Q,QT        = SEM_QQT(glnum)

lbc = true
rbc = false
GLSetBC!(L,lbc,rbc,ifperiodic)
GLSetBC!(AL,lbc,rbc,ifperiodic)

Lg    = QT*L*Q
ALg   = QT*AL*Q
Bg    = QT*B 



println("Done.")








# Testing the module

include("../Module_SEM1D/SEM1D.jl")
#using .SEM1D

#include("Module_StepperArnoldi/StepperArnoldi.jl")
#using .StepperArnoldi

#include("Module_CenterManifold/CenterManifold.jl")
#using .CenterManifold

using LinearAlgebra
using SparseArrays
using Printf
using PolynomialBases
using IterativeSolvers

using PyPlot

include("GL_Functions.jl")
#-------------------------------------------------- 

# Input parameters
Inp   = Get_SEM1D_Input()
# Nodal Bases
B0    = LobattoLegendre(Inp.N)
Bd    = LobattoLegendre(Inp.Nd)
# Geometric Matrices
GeoM  = SEM1D.SEMGeoMat(B0,Bd,Inp)

# GL Parameters
δ     = Set_GL_Params()
δc    = conj.(δ)

# GinzburgLandau Linear Operators
L,  B, OP,  Conv,  Src,  Lap  = SEM1D.GinzburgLandauSparse(δ, Inp,GeoM,B0)
LC, B, OPC, ConvC, SrcC, LapC = SEM1D.GinzburgLandauSparse(δc,Inp,GeoM,B0)

AL,  B, AOP, AConv,  ASrc,  ALap  = SEM1D.AdjointGinzburgLandauSparse(δ,Inp,GeoM,B0)
ALC, B,AOPC, AConvC, ASrcC, ALapC = SEM1D.AdjointGinzburgLandauSparse(δ,Inp,GeoM,B0)

ifperiodic  = false
ndof, glnum = SEM1D.SEM_Global_Num(GeoM.xm1,ifperiodic)
Q,QT        = SEM1D.SEM_QQT(glnum)
vmult       = Q*QT*ones(Float64,length(GeoM.xm1[:]))
vimult      = 1.0./vmult

SEM1D.GLSetBC!(L ,Inp.lbc,Inp.rbc,ifperiodic)
SEM1D.GLSetBC!(LC,Inp.lbc,Inp.rbc,ifperiodic)
SEM1D.GLSetBC!(AL,Inp.lbc,Inp.rbc,ifperiodic)
SEM1D.GLSetBC!(ALC,Inp.lbc,Inp.rbc,ifperiodic)

Lg    = QT*L*Q          # Linear Matrix
LCg   = QT*LC*Q         # Conjugated Linear Matrix

ALg   = QT*AL*Q         # Adjoint Linear Matrix
ALCg  = QT*ALC*Q        # Conjugated Adjoint Linear Matrix

Bg    = QT*B            # Mass Matrix
Bgi   = 1.0./Bg         # Inverse Mass Matrix
BgM   = diagm(Bg)
# OPg   = Bgi.*Lg         # Final Direct Operator
# OPCg  = Bgi.*LCg        # Final Conjugated Direct Operator
# AOPg  = Bgi.*ALg        # Final Adjoint Operator
# 
# Lapg  = Bgi.*(QT*(1.0/δ[4])*Lap*Q)        # Laplacian Operator
# LapCg = Bgi.*(QT*(1.0/δ[4]')*LapC*Q)      # Conjugated Laplacian Operator

OPg   = Lg         # Final Direct Operator
OPCg  = LCg        # Final Conjugated Direct Operator
AOPg  = ALg        # Final Adjoint Operator
AOPCg = ALCg       # Final Conjugated Adjoint Operator

Lapg  = (QT*(1.0/δ[4])*Lap*Q)        # Laplacian Operator
LapCg = (QT*(1.0/δ[4]')*LapC*Q)      # Conjugated Laplacian Operator

BiLapg  = diagm(Bgi)*(QT*(1.0/δ[4])*Lap*Q)        # (B^-1)*Laplacian Operator
BiLapCg = diagm(Bgi)*(QT*(1.0/δ[4]')*LapC*Q)      # (B^-1)*Conjugated Laplacian Operator

xg    = QT*(vimult.*GeoM.xm1[:])

println("Ginzburg Landau Setup Done.")
















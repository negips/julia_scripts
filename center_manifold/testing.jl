# Testing the module

include("Module_SEM1D/SEM1D.jl")
#using .SEM1D

include("Module_StepperArnoldi/StepperArnoldi.jl")
#using .StepperArnoldi

using LinearAlgebra
using SparseArrays
using Printf
using PolynomialBases


using PyPlot

N     = 12
Nd    = 19
nel   = 41
xs    = 0.0
xe    = 40.0
lbc   = true
rbc   = false

# Input parameters
Inp   = SEM1D.SEMInput(N,Nd,nel,xs,xe,lbc,rbc)
# Nodal Bases
B0    = LobattoLegendre(Inp.N)
Bd    = LobattoLegendre(Inp.Nd)
# Geometric Matrices
GeoM  = SEM1D.SEMGeoMat(B0,Bd,Inp)

δ     = ones(ComplexF64,4)    # Parameters
δ[1]  = -1.0                  # -U
δ[2]  =  0.741 + 1.025im      #  μ0
δ[3]  = -0.125                #  μx
δ[4]  = (1.0 - im)/sqrt(2.0)  #  γ
ifconj = false
if (ifconj)
  δ = con.(δ)
end  

# GinzburgLandau Linear Operators
L,  B, OP,  Conv,   Src, Lap  = SEM1D.GinzburgLandauSparse(δ,Inp,GeoM,B0)
AL, B, AOP, AConv, ASrc, ALap = SEM1D.AdjointGinzburgLandauSparse(δ,Inp,GeoM,B0)

ifperiodic = false
ndof, glnum = SEM1D.SEM_Global_Num(GeoM.xm1,ifperiodic)
Q,QT        = SEM1D.SEM_QQT(glnum)
vmult       = Q*QT*ones(Float64,length(GeoM.xm1[:]))
vimult      = 1.0./vmult

SEM1D.GLSetBC!(L ,Inp.lbc,Inp.rbc,ifperiodic)
SEM1D.GLSetBC!(AL,Inp.lbc,Inp.rbc,ifperiodic)

Lg    = QT*L*Q
ALg   = QT*AL*Q
Bg    = QT*B
Bgi   = 1.0./Bg
OPg   = Bgi.*Lg
AOPg  = Bgi.*ALg

# Stepper-Arnoldi
#-------------------------------------------------- 
ifadjoint         = false
ifoptimal         = false
ifverbose         = false
verbosestep       = 500
nsteps            = 500
dt                = 1.0e-4
StpInp            = StepperArnoldi.StepperInput(ifadjoint,ifoptimal,ifverbose,verbosestep,nsteps,dt)

ifarnoldi         = true 
ifverbose         = true
vlen              = ndof
nev               = 5
ekryl             = 15  
lkryl             = nev + ekryl 
ngs               = 2
bsize             = 1
outer_iterations  = 100
tol               = 1.0e-10
ArnInp            = StepperArnoldi.ArnoldiInput(ifarnoldi,ifverbose,vlen,nev,ekryl,lkryl,ngs,bsize,outer_iterations,tol)

Ω     = SEM1D.GLAnalyticalSpectra(δ)
if (StpInp.ifadjoint)
  Ω   = conj.(Ω)
end  

close("all")
h1    = figure(num=1,figsize=[8.,6.]);
ax1   = gca()
pl1   = ax1.plot(imag.(Ω),real.(Ω),linestyle="none",marker="o",markersize=8)
ax1.set_ylim([-2.75,0.5])
if (StpInp.ifadjoint)
  ax1.set_xlim([-1.80,0.0])
else  
  ax1.set_xlim([0.0,1.80])
end  

#ev2   = eigvals(Matrix(OPg))
#pl2   = ax1.plot(imag.(ev2),real.(ev2),linestyle="none",marker="x",markersize=4)

# Calculate Spectra
if StpInp.ifadjoint
  ArnOut      = StepperArnoldi.StepArn(AOPg,Bg,StpInp,ArnInp,Inp.lbc,Inp.rbc)
else  
  ArnOut      = StepperArnoldi.StepArn(OPg,Bg,StpInp,ArnInp,Inp.lbc,Inp.rbc)
end  

pl3   = ax1.plot(imag.(ArnOut.evals),real.(ArnOut.evals),linestyle="none",marker="s",markersize=4)

xg    = QT*(vimult.*GeoM.xm1[:])
h2    = figure(num=2,figsize=[8.,6.])
ax2   = gca()
for i in 1:ArnInp.nev
  ax2.plot(xg,real.(ArnOut.evecs[:,i]))
end  

println("Done.")








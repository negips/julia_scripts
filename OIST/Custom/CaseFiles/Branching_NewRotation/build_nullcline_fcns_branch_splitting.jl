#!/bin/julia

println("Building a null cline based on given Points/gradients by Minimizing the Lagrangian (Error)")

using LinearAlgebra
using Roots
using PyPlot
using Printf

const SRC = "/home/prabal/workstation/git/julia/OIST/Custom"

include("$SRC/BuildTimeDerivatives.jl")
include("$SRC/NullClines.jl")
include("$SRC/NullClineFcn.jl")
include("$SRC/ElementOf.jl")
include("$SRC/NullClineParams.jl")
include("$JULIACOMMON/MoveFigure.jl")

close("all")

lafs = 16

#include("select_nullclines.jl")
sets              = [207]     # 21

cm                = get_cmap("tab10")

set               = sets[1]
pars              = GetNullClineParams(set) 

h1                = figure(num=1)
ax1               = h1.subplots()
#MoveFigure(h1,10,10)
MoveFigure(h1,1250,830)

# Container to hold plot handles
PlotContainers    = Array{Any}(undef,20)


# Activator
#---------------------------------------------------------------------- 

f(x,y)            = FXY(x,y,pars.fc0,pars.fcx,pars.fcy)
if f(0.0,100.0)>0
  pars.fc0        = -pars.fc0
  pars.fcx        = -pars.fcx
  pars.fcy        = -pars.fcy
end  

# If we want +λ to be stabilizing or destabilizing
stabilizing = true
# Translated
if stabilizing 
  ϕfd   = 120 #150.0
else
  ϕfd   = -80.0
end  
println("F(x,y) Translated with Slope: $ϕfd Degrees")
ϕf      = ϕfd*π/180.0

# Plot the null-cline
λ0      = 0.0
dλ      = 0.6
#λvalues = [λ0-dλ; λ0; λ0+dλ]
λvalues = [λ0]
θvalues = [0.0]
plc = 0
for λ in λvalues
  local θ             = λ*pi/180.0
  # local f2(x,y)       = RotFXY(x,y,θ0,pars.fc0,pars.fcx,pars.fcy)
  local f2(x,y)       = TransFXY(x,y,λ,ϕf,pars.fc0,pars.fcx,pars.fcy)
  # Plot the null-cline
  local xi            = -20.0
  local yr0           = -10.0
  local yr1           =  50.0
  local dτ            = 1.0e-3
  local nsteps        = 200000
  local f0x,f0y       = NullClines(f2,xi,yr0,yr1,nsteps,dτ)
  global plc         += 1
  #PlotContainers[plc] = ax1.plot(f0x,f0y,linestyle="-",label="λ=$λ;")
  PlotContainers[plc] = ax1.plot(f0x,f0y,linestyle="-")
end  


# Inhibitor
#---------------------------------------------------------------------- 
g(x,y)            = FXY(x,y,pars.gc0,pars.gcx,pars.gcy)
if g(100.0,0.0)>0
  pars.gc0        = -pars.gc0
  pars.gcx        = -pars.gcx
  pars.gcy        = -pars.gcy
  g(x,y)          = FXY(x,y,pars.gc0,pars.gcx,pars.gcy)
end  
stabilizing = true
# Translated
if stabilizing
  ϕgd             = -30.0
else
  ϕgd             = 150.0
end  
#println("G(x,y) Translated with Slope: $ϕgd Degrees")
ϕg                = ϕgd*π/180.0

λvalues = [0.0]
θ0      =  0.0
dθ      =  -15.0/1.8
θvalues = [θ0; θ0+dθ*2; θ0-dθ*2]
Axis_X0 = 0.0
Axis_Y0 = -pars.gcx[1]/pars.gcy[1]*Axis_X0
#θvalues = [θ0]
for λ in θvalues
  local θ             = λ*pi/180.0
  # local g2(x,y)       = RotFXY(x,y,θ,pars.gc0,pars.gcx,pars.gcy)
  # local g2(x,y)       = RotXYFXY(x,y,Axis_X0,Axis_Y0,θ,pars.gc0,pars.gcx,pars.gcy)
  local g2(x,y)       = RotLinearFXY2(x,y,θ,pars.gc0,pars.gcx,pars.gcy)

  local xi            = -2.0
  local yr0           = -100.0
  local yr1           =  10.0
  local dτ            = 1.0e-3
  local nsteps        = 200000
  local g20x,g20y     = NullClines(g2,xi,yr0,yr1,nsteps,dτ)
  global plc         += 1

  if θ < 0.0
    PlotContainers[plc] = ax1.plot(g20x,g20y,linestyle="--",label="λ=2.0")
  elseif θ > 0.0 
    PlotContainers[plc] = ax1.plot(g20x,g20y,linestyle="--",label="λ=-2.0")
  else
    PlotContainers[plc] = ax1.plot(g20x,g20y,linestyle="-",label="λ=0.0")
  end

  # legend()
end  
legend()
# ax1.legend(PlotContainers[[2;3;4]],["λ=-2";"λ=0";"λ=2"])
#lg = ax1.legend([PlotContainers[2], PlotContainers[3], PlotContainers[4]], ["λ=-2", "λ=0", "λ=-2"]);
#legend()

#PlotContainers[6] = ax1.plot(pars.xB,pars.yB,linestyle=" ",marker="o",fillstyle="none")
#PlotContainers[7] = ax1.plot(pars.xdxB,pars.ydxB,linestyle=" ",marker="x")
#PlotContainers[8] = ax1.plot(pars.xdyB,pars.ydyB,linestyle=" ",marker="x")

ax1.set_xlabel(L"B", fontsize=lafs)
ax1.set_ylabel(L"A", fontsize=lafs)

ax1.set_xlim(-1.5,5.0)
ax1.set_ylim(-1.5,6.0)

MoveFigure(h1,1250,830)
fname0   = @sprintf "./plots/nullclines"
h1.savefig(fname0)

# Time dependent null-cline functions
yin     = LinRange(-1.5,7.0,5000)
#gt(z)   = -1.0/pars.gcx[1]*TransFXY(0.0,yin,z,ϕg,pars.gc0,pars.gcx,pars.gcy)
ft(z)   = -1.0/pars.fcx[1]*TransFXY(0.0,yin,z,ϕf,pars.fc0,pars.fcx,pars.fcy)

#gt(z)   = -1.0/pars.gcx[1]*RotFXY(0.0,yin,z,pars.gc0,pars.gcx,pars.gcy)
#gt(z)    = -1.0/pars.gcx[1]*RotXYFXY(0.0,yin,Axis_X0,Axis_Y0,z,pars.gc0,pars.gcx,pars.gcy)
slopeϕ  = atan(pars.gcx[1],pars.gcy[1])
gt(z)    = -1.0/(tan(slopeϕ-z)*pars.gcy[1])*RotLinearFXY2(0.0,yin,z,pars.gc0,pars.gcx,pars.gcy)
#ft(z)   = -1.0/pars.fcx[1]*RotFXY(0.0,yin,z,pars.fc0,pars.fcx,pars.fcy)

#PlotContainers[6] = ax1.plot(gt(-20.0*π/180.0),yin,linestyle="-.",linewidth=2,color=cm(5));

# Build Nullcline for the dynamic switching
#---------------------------------------- 
set               = 55
parsS             = GetNullClineParams(set)
δ                 = 0.005
λdot0(x,y)        = (1.0/δ)*FXY(x,y,parsS.fc0,parsS.fcx,parsS.fcy)
if λdot0(0.0,100.0)>0
  parsS.fc0        = -parsS.fc0
  parsS.fcx        = -parsS.fcx
  parsS.fcy        = -parsS.fcy
end  
λdot1(x,y)         = (1.0/δ)*FXY(x,y,parsS.fc0,parsS.fcx,parsS.fcy)

#α                 = 1.0
#xc                = 1.0
#λdot(x,y)         = (y-α*x)*(x^2 - xc^2)
xi                =  10.0
yr0               =  0.0
yr1               =  15.0
dτ                =  -1.0e-3
nsteps            = 50000
λdot0x1,λdot0y1   = NullClines(λdot1,xi,yr0,yr1,nsteps,dτ)

h4                = figure(num=4)
ax4               = h4.subplots()
ax4.plot(λdot0x1,λdot0y1,color=cm(3),linestyle="--")
ax4.set_ylabel(L"λ", fontsize=lafs)
ax4.set_xlabel(L"\widebar{A}", fontsize=lafs)

ax4.set_xlim(0.1,0.5)
ax4.set_ylim(-2.25,2.25)
fname0   = @sprintf "./plots/paramnullcline"
h4.savefig(fname0)

pause(0.01)

ϵ  = 0.1
η  = 1.0

println("Press x to stop. Any other key to continue")
xin = readline()
#xin = "x"
#xin = "y"
if xin !="x"
  # G(x,y,z) = TransFXY(x,y,z,ϕg,pars.gc0,pars.gcx,pars.gcy)*η
  F(x,y,z) = TransFXY(x,y,z,ϕf,pars.fc0,pars.fcx,pars.fcy)/ϵ

  # G(x,y,z)  = RotFXY(x,y,z,pars.gc0,pars.gcx,pars.gcy)*η
  #G(x,y,z)  = RotXYFXY(x,y,Axis_X0,Axis_Y0,z,pars.gc0,pars.gcx,pars.gcy)
  G(x,y,z)  = RotLinearFXY2(x,y,z,pars.gc0,pars.gcx,pars.gcy)
  # F(x,y,z)  = RotFXY(x,y,z,pars.fc0,pars.fcx,pars.fcy)/ϵ

  # Λ(x,y,z) = FXYZ(x,y,z,λc0,λcx,λcy,λcz)
  Flow(x,y,z1,z2) = [G(x,y,z1) F(x,y,z2)]

  close(h4)
  include("time_stepper_multiple_branch_splitting.jl")
end








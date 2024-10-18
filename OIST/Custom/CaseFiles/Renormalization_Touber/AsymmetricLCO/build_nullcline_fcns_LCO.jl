#!/bin/julia

println("Building a null cline based on given Points/gradients by Minimizing the Lagrangian (Error)")


using LinearAlgebra
using Roots
using PyPlot
using Printf

const SRC = "/home/prabal/workstation/git/julia/OIST/Custom"

include("../print_params.jl")
include("$SRC/BuildTimeDerivatives.jl")
include("$SRC/NullClines.jl")
include("$SRC/NullClineFcn.jl")
include("$SRC/ElementOf.jl")
include("$SRC/NullClineParams.jl")
include("$JULIACOMMON/MoveFigure.jl")

include("../renormalize_system.jl")

close("all")

lafs = 16

ifrenorm = true
iftouber = false

#include("select_nullclines.jl")
sets              = [204]

cm                = get_cmap("tab10")

# F,G
set               = sets[1]
pars              = GetNullClineParams(set) 
# λ
set               = 56
parsS             = GetNullClineParams(set)

if (ifrenorm) 
  Anorm,Bnorm     = renormalize_system!(pars)
  λnorm           = renormalize_λsystem!(parsS)
else
  Anorm           = 1.0
  Bnorm           = 1.0
  λnorm           = 1.0
end  

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
  local f2(x,y)       = TransFXY(x,y,λ,0.0,pars.fc0,pars.fcx,pars.fcy)
  # Plot the null-cline
  local xi            = -20.0
  local yr0           = -10.0
  local yr1           =  50.0
  local dτ            = 1.0e-3
  local nsteps        = 200000
  local f0x,f0y       = NullClines(f2,xi,yr0,yr1,nsteps,dτ)
  global plc         += 1
  PlotContainers[plc] = ax1.plot(f0x,f0y,linestyle="-",label="λ=$λ;")
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

λvalues = [0.0]
θ0      =  0.0
dθ      =  (-20.0)*λnorm
θvalues = [θ0]
Axis_X0 = 0.0
Axis_Y0 = -pars.gcx[1]/pars.gcy[1]*Axis_X0
#θvalues = [θ0]
for λ in θvalues
  local θ             = λ*pi/180.0
  #local g2(x,y)       = RotFXY(x,y,θ,pars.gc0,pars.gcx,pars.gcy)
  local g2(x,y)       = RotXYFXY(x,y,Axis_X0,Axis_Y0,θ,pars.gc0,pars.gcx,pars.gcy)
  local xi            = -2.0
  local yr0           = -100.0
  local yr1           =  10.0
  local dτ            = 1.0e-3
  local nsteps        = 200000
  local g20x,g20y     = NullClines(g2,xi,yr0,yr1,nsteps,dτ)
  global plc         += 1
  PlotContainers[plc] = ax1.plot(g20x,g20y,linestyle="--",label="θ=$λ")
  # legend()
end  
#legend()

if (ifrenorm)
  ax1.set_xlim(-2.0,8.0)
  ax1.set_ylim(-2.0,8.0)
else
  ax1.set_xlim(-1.5,5.0)
  ax1.set_ylim(-1.5,5.0)
end  

ax1.set_xlabel(L"B", fontsize=lafs)
ax1.set_ylabel(L"A", fontsize=lafs)

MoveFigure(h1,1250,830)
fname0   = @sprintf "./plots/nullclines"
h1.savefig(fname0)


# Time dependent null-cline functions
yin     = LinRange(-1.5,7.0,5000)
ft(z)   = -1.0/pars.fcx[1]*TransFXY(0.0,yin,z,0.0,pars.fc0,pars.fcx,pars.fcy)
gt(z)   = -1.0/pars.gcx[1]*RotXYFXY(0.0,yin,Axis_X0,Axis_Y0,z,pars.gc0,pars.gcx,pars.gcy)


# Build Nullcline for the dynamic switching
#---------------------------------------- 
δ                 = 0.0015
λdot0(x,y)        = (1.0/δ)*FXY(x,y,parsS.fc0,parsS.fcx,parsS.fcy)
if λdot0(0.0,100.0)>0
  parsS.fc0        = -parsS.fc0
  parsS.fcx        = -parsS.fcx
  parsS.fcy        = -parsS.fcy
end  
λdot1(x,y)         = (1.0/δ)*FXY(x,y,parsS.fc0,parsS.fcx,parsS.fcy)

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
if (ifrenorm) 
  ax4.set_xlim(-0.4,0.4)
  ax4.set_ylim(-1.6,1.6)
else  
  ax4.set_xlim(-0.4,0.4)
  ax4.set_ylim(-2.5,2.5)
end  

fname0   = @sprintf "./plots/paramnullcline"
h4.savefig(fname0)

pause(0.01)

ϵ  = 0.1
η  = 1.0
F(x,y,z)  = TransFXY(x,y,z,0.0,pars.fc0,pars.fcx,pars.fcy)/ϵ
G(x,y,z)  = RotXYFXY(x,y,Axis_X0,Axis_Y0,z,pars.gc0,pars.gcx,pars.gcy)

Flow(x,y,z1,z2) = [G(x,y,z1) F(x,y,z2)]


println("Press x to stop. Any other key to continue")
xin = readline()
# xin = "x"
#xin = "y"
if xin !="x"
  close(h4)
  include("time_stepper_multiple_LCO.jl")
end

print_params(F,G,λdot1,pars,parsS)

ndec = 5
@printf "ϵ   = %.*f\n" ndec ϵ
@printf "η   = %.*f\n" ndec η






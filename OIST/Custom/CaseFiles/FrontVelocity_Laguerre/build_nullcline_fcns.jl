#!/bin/julia

println("Building a null cline based on given Points/gradients by Minimizing the Lagrangian (Error)")


using LinearAlgebra
using Roots
using PyPlot

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
sets              = [2]

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

plc = 0
for i in 1:1
  # Plot the null-cline
  local xi            = -20.0
  local yr0           = -10.0
  local yr1           =  50.0
  local dτ            = 1.0e-3
  local nsteps        = 200000
  local f0x,f0y       = NullClines(f,xi,yr0,yr1,nsteps,dτ)
  global plc         += 1
  PlotContainers[plc] = ax1.plot(f0x,f0y,linestyle="-",label="f(x,y);")
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

for i in 1:1
  # Plot Null-cline
  local xi            = -10.0
  local yr0           = -30.0
  local yr1           =  10.0
  local dτ            = 1.0e-3
  local nsteps        = 200000
  local g20x,g20y     = NullClines(g,xi,yr0,yr1,nsteps,dτ)
  global plc         += 1
  PlotContainers[plc] = ax1.plot(g20x,g20y,linestyle="--",label="g(x,y)")
end  

ax1.set_xlabel(L"B", fontsize=lafs)
ax1.set_ylabel(L"A", fontsize=lafs)

ax1.set_xlim(-1.2,6.0)
ax1.set_ylim(-1.5,6.0)

MoveFigure(h1,1250,830)


# Time dependent null-cline functions
pause(0.01)

ϵ   = 0.1
η  = 1.0

println("Press x to stop. Any other key to continue")
#xin = readline()
#xin = "x"
xin = "y"
if xin !="x"
  F(x,y) = f(x,y)/ϵ
  G(x,y) = g(x,y)*η

  Flow(x,y) = [G(x,y) F(x,y)]

  include("sem_init_ref.jl")
  include("custom_params.jl")
  include("sem_main.jl")
  include("newton.jl")
#  include("time_stepper_multiple_oscillation.jl")
end








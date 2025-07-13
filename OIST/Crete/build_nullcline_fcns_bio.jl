#!/bin/julia

println("Building a null cline based on given Points/gradients by Minimizing the Lagrangian (Error)")


using LinearAlgebra
using Roots
using PyPlot
using Printf
using Random

const SRC = "/home/prabal/workstation/git/julia/OIST/Custom"

include("$SRC/BuildTimeDerivatives.jl")
include("$SRC/NullClines.jl")
include("$SRC/NullClineFcn.jl")
include("$SRC/ElementOf.jl")
include("$SRC/NullClineParams.jl")
include("$JULIACOMMON/MoveFigure.jl")

include("print_params.jl")
include("renormalize_system.jl")

close("all")

lafs = 16

ifrenorm = true
iftouber = false
ifλnorm2 = true

#include("select_nullclines.jl")
#sets              = [201]

cm                = get_cmap("tab10")

# Round off parameters
#round_params!(pars,parsS,4)

h1                = figure(num=1)
ax1               = h1.subplots()
#MoveFigure(h1,10,10)
MoveFigure(h1,1250,830)

# Container to hold plot handles
PlotContainers    = Array{Any}(undef,20)
plc               = 0

# dϕ/dt
gc0         = 0.0
gc1         = [-1.0]                 # ϕ
gc2         = [1.5]                  # X # 1.5
g(x,y)      = FXY(x,y,gc0,gc1,gc2)

# d²X/dt²
fc0         = 0.0
fc1         = [-5.0]                 # ϕ
fc2         = [-6.0; 7.0; -1.0]      # X
f(x,y)      = FXY(x,y,fc0,fc1,fc2)

xi          = -5.0
yr0         = -10.0
yr1         =  50.0
dτ          = 1.0e-3
nsteps      = 200000
g0x,g0y     = NullClines(g,xi,yr0,yr1,nsteps,dτ)
f0x,f0y     = NullClines(f,xi,yr0,yr1,nsteps,dτ)
plc        += 1
PlotContainers[plc] = ax1.plot(g0x,g0y,linestyle="-")
plc        += 1
PlotContainers[plc] = ax1.plot(f0x,f0y,linestyle="-")

ax1.set_xlim(-2.0,7.0)
ax1.set_ylim(-2.0,7.0)



MoveFigure(h1,1250,830)
fname0   = @sprintf "./plots/nullclines"
#h1.savefig(fname0)

ϵ  = 1.0E-1
η  = 1.0E+0
Flow(x,y)   = [g(x,y)/η f(x,y)/ϵ]

println("Press x to stop. Any other key to continue")
xin = readline()
#xin = "x"
#xin = "y"
if xin !="x"
  # close(h4)
  include("time_stepper_multiple_init.jl")
end










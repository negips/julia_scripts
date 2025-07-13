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
include("FitzhughNagumo.jl")
include("custom_params.jl")

include("print_params.jl")
include("renormalize_system.jl")

close("all")

lafs        = 16
cm          = get_cmap("tab10")
cmapU       = cm(0)
cmapV       = cm(1)

h1          = figure(num=1)
ax1         = h1.subplots()
MoveFigure(h1,1250,830)

# Container to hold plot handles
PlotContainers    = Array{Any}(undef,20)
plc               = 0

# Get the nullclines
Vx,Vy,Ux,Uy = FitzhughNagumoNullClines2(a,b,R,m)

plc        += 1
PlotContainers[plc] = ax1.plot(Vx,Vy,color=cmapV,linestyle="-")
plc        += 1
PlotContainers[plc] = ax1.plot(Ux,Uy,color=cmapU,linestyle="-")

ax1.set_xlim(-1.5,5.0)
ax1.set_ylim(-1.5,6.0)
ax1.set_xlabel(L"v", fontsize=lafs)
ax1.set_ylabel(L"u", fontsize=lafs)


MoveFigure(h1,1250,830)
fname0   = @sprintf "./plots/nullclines"
#h1.savefig(fname0)

# ϵ  = 1.0E-1
# η  = 1.0E+0

#FHN(x,y) = FitzhughNagumo2(x,y,a,b,c)
g(x,y)   = FHN2_G(x,y,R,m)
f(x,y)   = FHN2_F(x,y,a,b)

Flow(x,y)   = [g(x,y)/η f(x,y)/ϵ]

println("Press x to stop. Any other key to continue")
#xin = readline()
#xin = "x"
xin = "y"
if xin !="x"
  # close(h4)
  include("time_stepper_multiple_init.jl")
end










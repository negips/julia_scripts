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
sets              = [201]

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

# dv/dt
gc0         = 0.0
gc1         = [0.0]                  # v
gc2         = [-5.0]                 # x
gc3         = [-6.0; 7.0; -1.0]      # ϕ
g(x,y,z)    = FXYZ(x,y,z,gc0,gc1,gc2,gc3)

# dx/dt
fc0         = 0.0
fc1         = [1.0]                  # v
fc2         = [0.0]                  # x
fc3         = [0.0]                  # ϕ
f(x,y,z)    = FXYZ(x,y,z,fc0,fc1,fc2,fc3)

# dϕ/dt
hc0         = 0.0
hc1         = [0.0]                  # v
hc2         = [-1.0]                 # x
hc3         = [1.5]                  # ϕ
h(x,y,z)    = FXYZ(x,y,z,hc0,hc1,hc2,hc3)


xi          = -5.0
yr0         = -10.0
yr1         =  50.0
dτ          =  1.0e-3
nsteps      =  200000
g2(y,z)     = g(0.0,y,z)
g0x,g0y     = NullClines(g2,xi,yr0,yr1,nsteps,dτ)
h2(y,z)     = h(0.0,y,z)
h0x,h0y     = NullClines(h2,xi,yr0,yr1,nsteps,dτ)
plc        += 1
PlotContainers[plc] = ax1.plot(g0x,g0y,linestyle="-")
plc        += 1
PlotContainers[plc] = ax1.plot(h0x,h0y,linestyle="-")

ax1.set_xlim(-2.0,7.0)
ax1.set_ylim(-2.0,7.0)




MoveFigure(h1,1250,830)
fname0   = @sprintf "./plots/nullclines"
#h1.savefig(fname0)

ϵ  = 0.1
η  = 1.0
Flow(x,y,z)   = [g(x,y,z) f(x,y,z) h(x,y,z)]

println("Press x to stop. Any other key to continue")
# xin = readline()
xin = "x"
#xin = "y"
if xin !="x"
  # close(h4)
  # include("time_stepper_multiple_oscillation_pulse.jl")
end

# ndec = 4 
# fmt = Printf.Format("%s\t = %.$(ndec)e\n")
# Printf.format(stdout,fmt,"ϵ",ϵ) #
# Printf.format(stdout,fmt,"η",η) #
# Printf.format(stdout,fmt,"ν",δ) #









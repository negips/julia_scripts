#!/bin/julia

println("Building a null clines for Phase-Field")


using LinearAlgebra
using Roots
using PyPlot
using Printf

const SRC = "/home/prabal/workstation/git/julia/OIST/Custom"

include("$SRC/print_params.jl")
include("$SRC/NullClines.jl")
include("$SRC/NullClineFcn.jl")
include("$SRC/ElementOf.jl")
include("$SRC/NullClineParams.jl")
include("$JULIACOMMON/MoveFigure.jl")

#---------------------------------------------------------------------- 
"""
   function DFX(x,y)

     Fx = -(x + tanh(y))

"""
function DFX(x,y,ϵ)

  α  = 1.0/(tanh(1.0/ϵ))

  df = -α*tanh.(y/ϵ) .- x

  return df
end  
#---------------------------------------------------------------------- 
"""
   function DFY(x,y)

     Fy = -(y + tanh(x))

"""
function DFY(x,y,ϵ)

  α  = 1.0/(tanh(1.0/ϵ))

  df = -α*tanh.(x/ϵ) .- y

  return df
end  
#---------------------------------------------------------------------- 

close("all")

lafs  = 16
lgfs  = 12

cm                = get_cmap("tab10")

h1                = figure(num=1)
ax1               = h1.subplots()
MoveFigure(h1,1250,830)

# Container to hold plot handles
PlotContainers    = Array{Any}(undef,20)
plc = 0

# ϕ1
#---------------------------------------------------------------------- 
ϵ           = 0.05
α           = 1.0/(tanh(1.0/ϵ))
npts        = 1000000

ϕ2_tmp      = LinRange(-5.0,5.0,npts)
f_null      = -α*tanh.(ϕ2_tmp/ϵ)

plc         += 1
PlotContainers[plc] = ax1.plot(f_null,ϕ2_tmp,linestyle="-",label=L"f(ϕ_{1},ϕ_{2})=0")


# ϕ2
#---------------------------------------------------------------------- 
ϕ1_tmp      = LinRange(-5.0,5.0,npts)
g_null      = -α*tanh.(ϕ1_tmp/ϵ)

plc         += 1
PlotContainers[plc] = ax1.plot(ϕ1_tmp,g_null,linestyle="-",label=L"g(ϕ_{1},ϕ_{2})=0")


ax1.set_xlim(-2.0,2.0)
ax1.set_ylim(-2.0,2.0)

legend(fontsize=lgfs)
ax1.set_xlabel(L"ϕ_{1}",fontsize=lafs)
ax1.set_ylabel(L"ϕ_{2}",fontsize=lafs)


# MoveFigure(h1,1250,830)
# fname0   = @sprintf "./plots/nullclines"
# h1.savefig(fname0)
# 
# 
 
#ϵ  = 0.1
η  = 1.0
F(x,y)            = DFX(x,y,ϵ)/η
G(x,y)            = DFY(x,y,ϵ)/η
Flow(x,y)         = [F(x,y) G(x,y)]

println("Press x to stop. Any other key to continue")
#xin = readline()
xin = "x"
#xin = "y"
if xin !="x"
  #close(h4)
  include("time_stepper_multiple_phasefield.jl")
end








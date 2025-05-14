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
   function DFX(x,y,ϵ)

     Fx = -(x + α*tanh(y/ϵ))

"""
function DFX(x,y,ϵ)

  α  = 1.0/(tanh(1.0/ϵ))

  df = -α*tanh.(y/ϵ) .- x

  return df
end  
#---------------------------------------------------------------------- 
"""
   function DFY(x,y,ϵ)

     Fy = -(y + α*tanh(x/ϵ))

"""
function DFY(x,y,ϵ)

  α  = 1.0/(tanh(1.0/ϵ))

  df = -α*tanh.(x/ϵ) .- y

  return df
end  
#---------------------------------------------------------------------- 
"""
   function DFX2(x,y,ϵ)

     Fx = -0.5*((2x - 1) + α*tanh((2y-1)/ϵ))

"""
function DFX2(x,y,ϵ)

  α  = 1.0/(tanh(1.0/ϵ))

  df = -0.5*(α*tanh.((2.0*y .- 1.0)/ϵ) .+ 2.0*x .- 1.0)

  return df
end  
#---------------------------------------------------------------------- 
"""
   function DFY2(x,y,ϵ)

     Fy = -0.5*((2y - 1) + α*tanh((2x-1)/ϵ))

"""
function DFY2(x,y,ϵ)

  α  = 1.0/(tanh(1.0/ϵ))

  df = -0.5*(α*tanh.((2.0*x .- 1.0)/ϵ) .+ 2.0*y .- 1.0)

  return df
end  
#---------------------------------------------------------------------- 

ftype = 1

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

if (ftype == 1)
  f_null      = -α*tanh.(ϕ2_tmp/ϵ)
else
  f_null      = -0.5*(-1.0 .+ α*tanh.((2.0*ϕ2_tmp .- 1.0)/ϵ))
end  

plc         += 1
PlotContainers[plc] = ax1.plot(f_null,ϕ2_tmp,linestyle="-",label=L"f(ϕ_{1},ϕ_{2})=0")


# ϕ2
#---------------------------------------------------------------------- 
ϕ1_tmp      = LinRange(-5.0,5.0,npts)
if (ftype == 1)
  g_null      = -α*tanh.(ϕ1_tmp/ϵ)
else  
  g_null      = -0.5*(-1.0 .+ α*tanh.((2.0*ϕ1_tmp .- 1.0)/ϵ))
end  

plc         += 1
PlotContainers[plc] = ax1.plot(ϕ1_tmp,g_null,linestyle="-",label=L"g(ϕ_{1},ϕ_{2})=0")

if (ftype == 1)
  ax1.set_xlim(-2.0,2.0)
  ax1.set_ylim(-2.0,2.0)
else
  ax1.set_xlim(-0.5,1.5)
  ax1.set_ylim(-0.5,1.5)
end

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
if (ftype == 1)
  F(x,y)            = DFX(x,y,ϵ)/η
  G(x,y)            = DFY(x,y,ϵ)/η
else
  F(x,y)            = DFX2(x,y,ϵ)/η
  G(x,y)            = DFY2(x,y,ϵ)/η
end  
Flow(x,y)         = [F(x,y) G(x,y)]

println("Press x to stop. Any other key to continue")
xin = readline()
#xin = "x"
#xin = "y"
if xin !="x"
  #close(h4)
  include("time_stepper_multiple_phasefield.jl")
end








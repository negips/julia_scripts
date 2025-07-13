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
lgfs        = 12
txtfs       = 10
mksz        = 10
mkew        = 2.0
lw          = 2.0

cm          = get_cmap("tab10")
cmapU       = cm(0)
cmapV       = cm(1)

ifarrows          = true
ifannotate        = true
iffixedpt         = true
iftrajectory      = true

h1          = figure(num=1)
ax1         = h1.subplots()
MoveFigure(h1,1250,830)

# Container to hold plot handles
PlotContainers    = Array{Any}(undef,20)
plc               = 0

# Get the nullclines
Vx,Vy,Ux,Uy = FitzhughNagumoNullClines2(a,b,R,m)

plc        += 1
PlotContainers[plc] = ax1.plot(Vx,Vy,color=cmapV,linestyle="-",label="g(v,u)=0")
plc        += 1
PlotContainers[plc] = ax1.plot(Ux,Uy,color=cmapU,linestyle="-",label="f(v,u)=0")

ax1.set_xlim(-4.0,8.0)
ax1.set_ylim(-4.0,8.0)
ax1.set_xlabel(L"v", fontsize=lafs)
ax1.set_ylabel(L"u", fontsize=lafs)
legend(loc="best",fontsize=lgfs)


MoveFigure(h1,1250,830)
fname0   = @sprintf "./plots/nullclines0"
h1.savefig(fname0)

if (ifarrows)
  # u-flow
  arr1 = arrow(-1.0,  7.5, 0.0,-1.0,width=0.1,color=cmapU)
  fname1   = @sprintf "./plots/nullclines1"
  h1.savefig(fname1)
  
  arr2 = arrow( 3.0, -3.5, 0.0, 1.0,width=0.1,color=cmapU)
  fname2   = @sprintf "./plots/nullclines2"
  h1.savefig(fname2)
  
  arr3 = arrow( 1.5, 2.75, 0.0, 0.5,width=0.1,color=cmapU)
  fname3   = @sprintf "./plots/nullclines3"
  h1.savefig(fname3)
  
  arr4 = arrow( 1.5, 1.5, 0.0, -0.5,width=0.1,color=cmapU)
  fname4   = @sprintf "./plots/nullclines4"
  h1.savefig(fname4)
  
  # v-flow
  arr5 = arrow(-1.0,  3.0, 1.0, 0.0,width=0.1,color=cmapV)
  fname5   = @sprintf "./plots/nullclines5"
  h1.savefig(fname5)
  
  arr6 = arrow( 5.0,  1.0,-1.0, 0.0,width=0.1,color=cmapV)
  fname6   = @sprintf "./plots/nullclines6"
  h1.savefig(fname6)
end  

if ifannotate
  
  txt1 = ax1.text(-3.0, 5.80, "Upper Branch (Attracting)", fontsize=txtfs,
                  rotation=-7, rotation_mode="anchor")
  txt2 = ax1.text(2.0, -1.30, "Lower Branch (Attracting)", fontsize=txtfs,
                  rotation=-7, rotation_mode="anchor")
  txt3 = ax1.text(-0.0,1.5, "Repelling", fontsize=txtfs,
                  rotation=25, rotation_mode="anchor") 
  txt4 = ax1.text(-3.0,-1.9, "Attracting", fontsize=txtfs,
                  rotation=30, rotation_mode="anchor")
  # Phase plot 
  annfname   = @sprintf "./plots/branches"
  h1.savefig(annfname)
end  

if iffixedpt
  plc        += 1
  PlotContainers[plc] = ax1.plot(0.0,0.0,color="black",linestyle="none",marker="o",markersize=mksz,fillstyle="none",markeredgewidth=mkew,label="Fixed Point.")
  fpfname   = @sprintf "./plots/fixedpoint"
  h1.savefig(fpfname)
end

# Trajectory
if (iftrajectory)
  g(x,y)   = FHN2_G(x,y,R,m)
  f(x,y)   = FHN2_F(x,y,a,b)

  trcol = "gray"
  include("phase_flow_fitzhughnagumo.jl")

  trfname1 = @sprintf "./plots/trajectory1"
  h1.savefig(trfname1)

  g(x,y)   = FHN2_G(x,y,R,m)/η
  f(x,y)   = FHN2_F(x,y,a,b)/ϵ

  trcol = "black"
  include("phase_flow_fitzhughnagumo.jl")

  trfname2 = @sprintf "./plots/trajectory2"
  h1.savefig(trfname2)

end










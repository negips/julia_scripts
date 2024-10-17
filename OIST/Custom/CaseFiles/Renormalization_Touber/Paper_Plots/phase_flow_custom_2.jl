#!/bin/julia

println("Building a null cline based on given Points/gradients by Minimizing the Lagrangian (Error)")


using LinearAlgebra
using Roots
using PyPlot
using Printf

const SRC = "/home/prabal/workstation/git/julia/OIST/Custom"

include("$SRC/print_params.jl")
# include("$SRC/GetDynamicNullCline.jl")
# include("$SRC/BuildTimeDerivatives.jl")
include("$SRC/NullClines.jl")
include("$SRC/NullClineFcn.jl")
include("$SRC/ElementOf.jl")
include("$SRC/NullClineParams.jl")
include("$JULIACOMMON/MoveFigure.jl")

include("../renormalize_system.jl")


function FXY(x,y,c0,cx,cy)
  nx = length(cx)
  ny = length(cy)

  s = c0
  for i in 1:nx
    s = s .+ cx[i]*x.^i
  end

  for j in 1:ny
    s = s .+ cy[j]*y.^j
  end  

  return s
end  
#---------------------------------------------------------------------- 
function FX(x,c0,cx)
  nx = length(cx)

  s = c0
  for i in 1:nx
    s = s .+ cx[i]*x.^i
  end

  return s
end  
#---------------------------------------------------------------------- 
function FY(y,c0,cy)
  ny = length(cy)

  s = c0
  for i in 1:ny
    s = s .+ cy[i]*y.^i
  end

  return s
end  
#---------------------------------------------------------------------- 

close("all")
lafs              = 16

ifrenorm = true
iftouber = false

#include("select_nullclines.jl")
sets              = [201]

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

xi                = -20.0
yr0               = -10.0
yr1               =  50.0
dτ                = 1.0e-3
nsteps            = 200000
f(x,y)            = FXY(x,y,pars.fc0,pars.fcx,pars.fcy)
if f(0.0,100.0)>0
  pars.fc0        = -pars.fc0
  pars.fcx        = -pars.fcx
  pars.fcy        = -pars.fcy
  f(x,y)          = FXY(x,y,pars.fc0,pars.fcx,pars.fcy)
end  
f0x,f0y           = NullClines(f,xi,yr0,yr1,nsteps,dτ)

# Container to hold plot handles
PlotContainers    = Array{Any}(undef,10)

PlotContainers[1] = ax1.plot(f0x,f0y,color=cm(0),label="f(A,B)=0")
# PlotContainers[2] = ax1.plot(pars.xA,pars.yA,linestyle=" ",marker="s",fillstyle="none")
# PlotContainers[3] = ax1.plot(pars.xdxA,pars.ydxA,linestyle=" ",marker="x")
# PlotContainers[4] = ax1.plot(pars.xdyA,pars.ydyA,linestyle=" ",marker="x")


xi                = -10.0
yr0               = -50.0
yr1               =  00.0
dτ                =  1.0e-3
nsteps            = 100000
g(x,y)            = FXY(x,y,pars.gc0,pars.gcx,pars.gcy)
if g(100.0,0.0)>0
  pars.gc0        = -pars.gc0
  pars.gcx        = -pars.gcx
  pars.gcy        = -pars.gcy
  g(x,y)          = FXY(x,y,pars.gc0,pars.gcx,pars.gcy)
end  

g0x,g0y           = NullClines(g,xi,yr0,yr1,nsteps,dτ)


PlotContainers[5] = ax1.plot(g0x,g0y,color=cm(1),label="g(A,B)=0")
# PlotContainers[6] = ax1.plot(pars.xB,pars.yB,linestyle=" ",marker="o",fillstyle="none")
# PlotContainers[7] = ax1.plot(pars.xdxB,pars.ydxB,linestyle=" ",marker="x")
# PlotContainers[8] = ax1.plot(pars.xdyB,pars.ydyB,linestyle=" ",marker="x")

if (ifrenorm)
  ax1.set_xlim(-3.5,10.0)
  ax1.set_ylim(-3.5,10.0)
else
  ax1.set_xlim(-1.5,5.0)
  ax1.set_ylim(-1.5,5.0)
end  
ax1.set_xlabel(L"B", fontsize=lafs)
ax1.set_ylabel(L"A", fontsize=lafs)

lg = legend(loc="center right")

#ax1.set_xlim(-3.0,8.0)
#ax1.set_ylim(-3.0,8.0)

MoveFigure(h1,1250,500)

pause(0.01)

println("Press x to stop. Any other key to continue")
xin = readline()
#xin = "x"
#xin = " "
ϵ  = 0.1
η  = 1.0
if xin !="x"
  F(x,y) =  f(x,y)/ϵ
  G(x,y) =  g(x,y)*η

  Flow(x,y) = [G(x,y) F(x,y)]
  include("phase_flow.jl")
end



#!/bin/julia

println("Building a null cline based on given Points/gradients by Minimizing the Lagrangian (Error)")


using LinearAlgebra
using Roots
using PyPlot

include("NullClines.jl")
include("NullClineFcn.jl")
include("ElementOf.jl")
include("NullClineParams.jl")
include("MoveFigure.jl")

#---------------------------------------------------------------------- 
function FXYZ(x,y,z,c0,cx,cy,cz)
  nx = length(cx)
  ny = length(cy)
  nz = length(cz)

  s = c0
  for i in 1:nx
    s = s .+ cx[i]*x.^i
  end

  for j in 1:ny
    s = s .+ cy[j]*y.^j
  end  

  for k in 1:nz
    s = s .+ cz[k]*z.^k
  end  
 
  return s
end  
#---------------------------------------------------------------------- 

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
function TransFXY(x,y,λ,θ,c0,cx,cy)

  nx = length(cx)
  ny = length(cy)

  x1 = x .- λ*cos(θ)    
  y1 = y .- λ*sin(θ)

  s = c0
  for i in 1:nx
    s = s .+ cx[i]*x1.^i
  end

  for j in 1:ny
    s = s .+ cy[j]*y1.^j
  end  

  return s
end  
#---------------------------------------------------------------------- 

close("all")

lafs = 16

include("select_nullclines.jl")

cm                = get_cmap("tab10")

set               = sets[1]
pars              = GetNullClineParams(set) 

h1                = figure(num=1)
ax1               = h1.subplots()
#MoveFigure(h1,10,10)
MoveFigure(h1,1250,830)

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

PlotContainers[1] = ax1.plot(f0x,f0y,color=cm(0))
PlotContainers[2] = ax1.plot(pars.xA,pars.yA,linestyle=" ",marker="s",fillstyle="none")
PlotContainers[3] = ax1.plot(pars.xdxA,pars.ydxA,linestyle=" ",marker="x")
PlotContainers[4] = ax1.plot(pars.xdyA,pars.ydyA,linestyle=" ",marker="x")

xi                = -10.0
yr0               = -20.0
yr1               =  10.0
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


PlotContainers[5] = ax1.plot(g0x,g0y,color=cm(1))
PlotContainers[6] = ax1.plot(pars.xB,pars.yB,linestyle=" ",marker="o",fillstyle="none")
PlotContainers[7] = ax1.plot(pars.xdxB,pars.ydxB,linestyle=" ",marker="x")
PlotContainers[8] = ax1.plot(pars.xdyB,pars.ydyB,linestyle=" ",marker="x")

ax1.set_xlabel(L"B", fontsize=lafs)
ax1.set_ylabel(L"A", fontsize=lafs)

ax1.set_xlim(-1.2,6.0)
ax1.set_ylim(-1.5,6.0)


# If we want λ to be stabilizing or destabilizing
stabilizing = true
# Translated
if stabilizing 
  ϕfd               = 150 #150.0
else
  ϕfd               = -30.0
end  

println("F(x,y) Translated with Slope: $ϕfd Degrees")
ϕf                = ϕfd*π/180.0
λ0                =  1.0
f2(x,y)           = TransFXY(x,y,λ0,ϕf,pars.fc0,pars.fcx,pars.fcy)

xi                = -20.0
yr0               = -10.0
yr1               =  50.0
dτ                = 1.0e-3
nsteps            = 200000
f20x,f20y         = NullClines(f2,xi,yr0,yr1,nsteps,dτ)

PlotContainers[9] = ax1.plot(f20x,f20y,color=cm(0),linestyle="--",label="λ=$λ0; ϕ=$ϕfd")


# Translated
if stabilizing
  ϕgd               = -30.0
else
  ϕgd               = 150.0
end  

println("G(x,y) Translated with Slope: $ϕgd Degrees")
ϕg                = ϕgd*π/180.0
#λ                 = 1.0
g2(x,y)           = TransFXY(x,y,λ0,ϕg,pars.gc0,pars.gcx,pars.gcy)

xi                = -20.0
yr0               = -30.0
yr1               =  10.0
dτ                = 1.0e-3
nsteps            = 200000
g20x,g20y         = NullClines(g2,xi,yr0,yr1,nsteps,dτ)

PlotContainers[10] = ax1.plot(g20x,g20y,color=cm(1),linestyle="--",label="λ=$λ0; ϕ=$ϕgd")
#legend()

#MoveFigure(h1,10,10)
MoveFigure(h1,1250,830)


# Time dependent null-cline functions
yin     = LinRange(-1.5,6.0,5000)
gt(z)   = -1.0/pars.gcx[1]*TransFXY(0.0,yin,z,ϕg,pars.gc0,pars.gcx,pars.gcy)
ft(z)   = -1.0/pars.fcx[1]*TransFXY(0.0,yin,z,ϕf,pars.fc0,pars.fcx,pars.fcy)


pause(0.01)

ϵ   = 0.1
η  = 1.0

# Equilibrium values of λ and a
eq_λ   = 10.0
eq_a   = 5.0
λc0    =  0.0                 # constant term
λcx    =  zeros(Float64,3)    # dependence on inhibitor b
λcy    =  zeros(Float64,3)    # dependence on activator a
λcz    =  zeros(Float64,3)    # dependence on λ
λc0    = 0.0
λcy[1] = -0.0/eq_a
λcy[2] = (1.0/eq_a)^2
λcz[1] = -(1.0/eq_λ)

println("Press x to stop. Any other key to continue")
xin = readline()
#xin = "x"
if xin !="x"
  G(x,y,z) = TransFXY(x,y,z,ϕg,pars.gc0,pars.gcx,pars.gcy)*η
  F(x,y,z) = TransFXY(x,y,z,ϕf,pars.fc0,pars.fcx,pars.fcy)/ϵ
  Λ(x,y,z) = FXYZ(x,y,z,λc0,λcx,λcy,λcz)
  Flow(x,y,z) = [G(x,y,z) F(x,y,z) Λ(x,y,z)]

#  include("time_stepper_multiple_init.jl")

  include("time_stepper_multiple_three.jl")
end








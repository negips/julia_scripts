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

  x1 = x - λ*cos(θ)    
  y1 = y - λ*sin(θ)

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

# Sets: 1         : Pulses
#     : 2         : Slugs. No change in instability threshold
#     : 3         : Pulses. Smaller instability threshold
#     : 4         : Slugs. Smaller instability threshold
#     : 5         : Extreme Slugs. Smaller instability threshold
#     : 6         : Extreme Pulse collapse
#     : 10        : Limit-cycline Oscillation. Activation dominated
#     : 11        : Extreme slugs
#     : 12        : Limit-cycline Oscillation. De-activation dominated
#     : 13        : Two fixed points - upper and lower branch.
#     : 14        : Symmetric fixed points - upper and lower branch
#     : 15        : Symmetric LCO
#     : 16        : Two Asymmetric fixed points
#     : 17        : Unstable G-null-cline
#     : 21        : Dynamic Switching (λ)
#     : 22        : Dynamic Switching (λ)

set               = 1
pars              = GetNullClineParams(set) 

h1                = figure(num=1)
ax1               = h1.subplots()

xi                = -20.0
yr0               = -10.0
yr1               =  50.0
dτ                = 1.0e-3
nsteps            = 200000
f(x,y)            = FXY(x,y,pars.fc0,pars.fcx,pars.fcy)
f0x,f0y           = NullClines(f,xi,yr0,yr1,nsteps,dτ)

ax1.plot(f0x,f0y,color=cm(0))
ax1.plot(pars.xA,pars.yA,linestyle=" ",marker="s",fillstyle="none")
ax1.plot(pars.xdxA,pars.ydxA,linestyle=" ",marker="x")
ax1.plot(pars.xdyA,pars.ydyA,linestyle=" ",marker="x")

xi                = -10.0
yr0               = -20.0
yr1               =  50.0
dτ                =  1.0e-3
nsteps            = 100000
g(x,y)            = FXY(x,y,pars.gc0,pars.gcx,pars.gcy)
g0x,g0y           = NullClines(g,xi,yr0,yr1,nsteps,dτ)


ax1.plot(g0x,g0y,color=cm(1))
ax1.plot(pars.xB,pars.yB,linestyle=" ",marker="o",fillstyle="none")
ax1.plot(pars.xdxB,pars.ydxB,linestyle=" ",marker="x")
ax1.plot(pars.xdyB,pars.ydyB,linestyle=" ",marker="x")

ax1.set_xlabel(L"B", fontsize=lafs)
ax1.set_ylabel(L"A", fontsize=lafs)

ax1.set_xlim(-1.2,6.0)
ax1.set_ylim(-1.5,6.0)

# Translated
ϕfd               = 150.0
println("F(x,y) Translated with Slope: $ϕfd Degrees")
ϕf                = ϕfd*π/180.0
λ                 = 1.0
f2(x,y)           = TransFXY(x,y,λ,ϕf,pars.fc0,pars.fcx,pars.fcy)

xi                = -20.0
yr0               = -10.0
yr1               =  50.0
dτ                = 1.0e-3
nsteps            = 200000
f20x,f20y           = NullClines(f2,xi,yr0,yr1,nsteps,dτ)

ax1.plot(f20x,f20y,color=cm(0),linestyle="--",label="λ=$λ; ϕ=$ϕfd")


# Translated
ϕgd               = 30.0
println("G(x,y) Translated with Slope: $ϕgd Degrees")
ϕg                = ϕg*π/180.0
λ                 = 1.0
g2(x,y)           = TransFXY(x,y,λ,ϕg,pars.gc0,pars.gcx,pars.gcy)

xi                = -20.0
yr0               = -30.0
yr1               =  10.0
dτ                = 1.0e-3
nsteps            = 200000
g20x,g20y         = NullClines(g2,xi,yr0,yr1,nsteps,dτ)

ax1.plot(g20x,g20y,color=cm(1),linestyle="--",label="λ=$λ; ϕ=$ϕgd")
legend()

MoveFigure(h1,10,10)

pause(0.01)

println("Press x to stop. Any other key to continue")
#xin = readline()
xin = "x"
if xin !="x"
  ϵ           = 1.0
  if f(0.0,100.0)>0
    F(x,y) = -f(x,y)/ϵ
  else
    F(x,y) =  f(x,y)/ϵ
  end  

  η  = 1.0
  if g(100.0,0.0) > 0
    G(x,y) = -η*g(x,y)
  else
    G(x,y) =  η*g(x,y)
  end
  Flow(x,y) = [G(x,y) F(x,y)]
  include("time_stepper_multiple.jl")
end


#!/bin/julia

println("Building a null cline based on given Points/gradients by Minimizing the Lagrangian (Error)")


using LinearAlgebra
using Roots
using PyPlot

include("NullClines.jl")
include("NullClineFcn.jl")
include("ElementOf.jl")
include("NullClineParams.jl")

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
function TERR_F(x,y,c0,cx,cy)

  n  = length(x)
  nx = length(cx)
  ny = length(cy)

  er = 0.0
  for k in 1:n
    s = c0
    for i in 1:nx
      s = s + cx[i]*x[k]^i
    end

    for j in 1:ny
      s = s + cy[j]*y[k]^j
    end
    er = er + s^2
  end  

  return er
end  
#---------------------------------------------------------------------- 

#ion()       # Show plots immediately

close("all")

lafs = 16

# Sets: 1         : Forward-Backward Traveling pulses
#     : 2         : Forward-Backward Traveling Slugs
#     : 3         : Limit-cycline Oscillation
#     : 4         : Testing. Trying to figure out Crossings.    


set               = 1
pars              = GetNullClineParams(set) 

h1                = figure(num=1)
ax1               = h1.subplots()

xi                = -20.0
yr0               = -50.0
yr1               = 50.0
dτ                = 1.0e-3
nsteps            = 120000
f(x,y)            = FXY(x,y,pars.fc0,pars.fcx,pars.fcy)
f0x,f0y           = NullClines(f,xi,yr0,yr1,nsteps,dτ)

ax1.plot(f0x,f0y)
ax1.plot(pars.xA,pars.yA,linestyle=" ",marker="s",fillstyle="none")
ax1.plot(pars.xdxA,pars.ydxA,linestyle=" ",marker="x")
ax1.plot(pars.xdyA,pars.ydyA,linestyle=" ",marker="x")

xi  = -10.0
yr0 = -50.0
yr1 = 1.0
dτ  = 1.0e-3
nsteps = 100000
g(x,y)            = FXY(x,y,pars.gc0,pars.gcx,pars.gcy)
g0x,g0y           = NullClines(g,xi,yr0,yr1,nsteps,dτ)




ax1.plot(g0x,g0y)
ax1.plot(pars.xB,pars.yB,linestyle=" ",marker="o",fillstyle="none")
ax1.plot(pars.xdxB,pars.ydxB,linestyle=" ",marker="x")
ax1.plot(pars.xdyB,pars.ydyB,linestyle=" ",marker="x")

ax1.set_xlabel(L"x", fontsize=lafs)
ax1.set_ylabel(L"y", fontsize=lafs)


pause(0.01)

println("Press x to stop. Any other key to continue")
xin = readline()
#xin = "x"
if xin !="x"
  ϵ           = 0.1
  if f(0.0,100.0)>0
    F(x,y) = -f(x,y)/ϵ
  else
    F(x,y) =  f(x,y)/ϵ
  end  

  η  = 1.0
  if g(100.0,0.0)>0
    G(x,y) = -η*g(x,y)
  else
    G(x,y) =  η*g(x,y)
  end
  Flow(x,y) = [G(x,y) F(x,y)]
  include("time_stepper_multiple.jl")
end




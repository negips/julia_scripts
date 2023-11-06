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
function RotFXY(x,y,Rc,Rs,c0,cx,cy)
  nx = length(cx)
  ny = length(cy)

  s = c0
  for i in 1:nx
    s = s .+ cx[i]*(Rc*x .+ Rs*y).^i
  end

  for j in 1:ny
    s = s .+ cy[j]*(-Rs*x .+ Rc*y).^j
  end  

  return s
end  
#---------------------------------------------------------------------- 
function RotXYFXY(x,y,Rc,Rs,x0,y0,c0,cx,cy)
  nx = length(cx)
  ny = length(cy)

  s = c0
  for i in 1:nx
    s = s .+ cx[i]*(x0 + Rc*(x-x0) .+ Rs*(y-y0)).^i
  end

  for j in 1:ny
    s = s .+ cy[j]*(y0 + -Rs*(x-x0) .+ Rc*(y-y0)).^j
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
<<<<<<< HEAD
#     : 15        : Symmetric LCO
#     : 16        : Two Asymmetric fixed points
#     : 21        : Dynamic Switching (λ)
#     : 22        : Dynamic Switching (λ)
#
#
sets              = [3 1]
=======
sets              = [5 3]
>>>>>>> master
nsets             = length(sets)
h1                = figure(num=1)
ax1               = h1.subplots()


set               = sets[1]
pars1             = GetNullClineParams(set) 
f1(x,y)           = FXY(x,y,pars1.fc0,pars1.fcx,pars1.fcy)
g1(x,y)           = FXY(x,y,pars1.gc0,pars1.gcx,pars1.gcy)

set               = sets[2]
pars2             = GetNullClineParams(set) 
f2(x,y)           = FXY(x,y,pars2.fc0,pars2.fcx,pars2.fcy)
g2(x,y)           = FXY(x,y,pars2.gc0,pars2.gcx,pars2.gcy)

xi                = -10.0
yr0               = -50.0
yr1               = 50.0
dτ                = 1.0e-3
nsteps            = 120000

f0x1,f0y1         = NullClines(f1,xi,yr0,yr1,nsteps,dτ)
f0x2,f0y2         = NullClines(f2,xi,yr0,yr1,nsteps,dτ)


cm    = get_cmap("tab10");

ax1.plot(f0x1,f0y1,color=cm(0))
ax1.plot(f0x2,f0y2,linestyle="--",color=cm(0))

xi  = -10.0
yr0 = -50.0
yr1 = 1.0
dτ  = 1.0e-3
nsteps = 100000
g0x1,g0y1          = NullClines(g1,xi,yr0,yr1,nsteps,dτ)
g0x2,g0y2          = NullClines(g2,xi,yr0,yr1,nsteps,dτ)

ax1.plot(g0x1,g0y1,color=cm(1))
ax1.plot(g0x2,g0y2,linestyle="--",color=cm(1))

ax1.set_xlabel(L"B", fontsize=lafs)
ax1.set_ylabel(L"A", fontsize=lafs)

#ax1.set_xlim(-1.2,3.0)
#ax1.set_ylim(-0.28,1.3)
ax1.set_xlim(-1.2,6.0)
ax1.set_ylim(-1.5,6.0)
MoveFigure(h1,1250,500)

# Build Nullcline for the dynamic switching
#---------------------------------------- 
set               = 23
parsS             = GetNullClineParams(set) 
λdot1(x,y)        = FXY(x,y,parsS.fc0,parsS.fcx,parsS.fcy)
#α                 = 1.0
#xc                = 1.0
#λdot(x,y)         = (y-α*x)*(x^2 - xc^2)
xi                = -10.0
yr0               =  0.0
yr1               =  15.0
dτ                =  1.0e-3
nsteps            = 100000
λdot0x1,λdot0y1   = NullClines(λdot1,xi,yr0,yr1,nsteps,dτ)

h4                = figure(num=4)
ax4               = h4.subplots()
#ax4.plot(λdot0x1,λdot0y1,color=cm(3),linestyle="--")
#ax4.set_xlim(-2.0,4.0)
#ax4.set_ylim(-2.0,2.0)
#ax4.grid()


ϕ                 = 1.3
println("Rotation Angle: $ϕ Degrees")
ϕ                 = ϕ*π/180.0
Rc                = cos(ϕ)
Rs                = sin(ϕ)
#λdot2(x,y)        = RotFXY(x,y,Rc,Rs,parsS.fc0,parsS.fcx,parsS.fcy)
λdot2(x,y)        = RotXYFXY(x,y,Rc,Rs,1.60,0.21,parsS.fc0,parsS.fcx,parsS.fcy)
λdot3(x,y)        = λdot2(x,y)
xi                = -10.0
yr0               = -5.0
yr1               = 15.0
dτ                =  1.0e-3
nsteps            = 100000
λdot0x1,λdot0y1   = NullClines(λdot3,xi,yr0,yr1,nsteps,dτ)
ax4.plot(λdot0x1,λdot0y1,color=cm(3),linestyle="-")
ax4.set_xlabel(L"θ", fontsize=lafs)
ax4.set_ylabel(L"λ", fontsize=lafs)

#ax4.set_xlim(0.0,20.0)
#ax4.set_ylim(-0.1,1.25)
ax4.grid()
MoveFigure(h4,600,500)

<<<<<<< HEAD
ϵλ    = 0.1
tmp1  = λdot0x1 .>= - 0.0 .&& λdot0x1 .< 7.0
=======
ϵλ    = 0.01
tmp1  = λdot0x1 .> - 2.0 .&& λdot0x1 .< 7.0
>>>>>>> master
tmp2  = λdot0y1.*tmp1
λmax  = maximum(tmp2)
if (λdot3(0.0,10.0)>0)
  λdot(x,y) = -λdot3(x,y)/ϵλ/λmax
else
  λdot(x,y) = λdot3(x,y)/ϵλ/λmax
end  

pause(0.01)

println("Press x to stop. Any other key to continue")
xin = readline()
#xin = "x"
if xin !="x"
  ϵ           = 0.1
  if f1(0.0,100.0)>0
    F1(x,y) = -f1(x,y)/ϵ
  else
    F1(x,y) =  f1(x,y)/ϵ
  end  
  if f2(0.0,100.0)>0
    F2(x,y) = -f2(x,y)/ϵ
  else
    F2(x,y) =  f2(x,y)/ϵ
  end  

  η  = 1.0
  if g1(100.0,0.0)>0
    G1(x,y) = -η*g1(x,y)
  else
    G1(x,y) =  η*g1(x,y)
  end
  if g2(100.0,0.0)>0
    G2(x,y) = -η*g2(x,y)
  else
    G2(x,y) =  η*g2(x,y)
  end

  Flow1(x,y) = [G1(x,y) F1(x,y)]
  Flow2(x,y) = [G2(x,y) F2(x,y)]

#  include("time_stepper_multiple_two.jl")
#  include("time_stepper_multiple_dynamic.jl")
  include("time_stepper_multiple_dynamic2.jl")
end




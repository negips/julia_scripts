#!/bin/julia

println("Building a null cline based on given Points/gradients by Minimizing the Lagrangian (Error)")


using LinearAlgebra
using Roots
using PyPlot

include("NullClines.jl")
include("NullClineFcn.jl")
include("ElementOf.jl")

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

nx = 1
ny = 3
n  = nx + ny + 1      # No of free parameters/Conditions to satisfy

FAC = 0.2

m = 3       # No of incident points
x = zeros(Float64,m)
y = zeros(Float64,m)
if m > 0
  x     = FAC*[0.0 -5.0 -5.0]
  y     = FAC*[0.0  1.0  6.0]
end  

# X-Derivative Points
mx = 0
Xx = zeros(Float64,mx)
Yx = zeros(Float64,mx)
if (mx>0)
  Xx = x[2]
  Yx = y[2] 
end  

# Y-Derivative Points
my = 1
Xy = zeros(Float64,my)
Yy = zeros(Float64,my)
if my>0
  Xy[1] = x[2]
  Yy[1] = y[2]
end  

fc = NullClineFcn(x,y,Xx,Yx,Xy,Yy,nx,ny)

#c = M\rhs
fc0 = fc[1]
fcx = fc[2:nx+1]
fcy = fc[nx+2:nx+ny+1]

f(x,y) = FXY(x,y,fc0,fcx,fcy)

xi = -20.0
yr0 = -50.0
yr1 = 50.0
dτ  = 1.0e-3
nsteps = 120000

xx,yy = NullClines(f,xi,yr0,yr1,nsteps,dτ)

h1  = figure(num=1)
ax1 = h1.subplots()

ax1.plot(xx,yy)
ax1.set_xlabel(L"x", fontsize=lafs)
ax1.set_ylabel(L"y", fontsize=lafs)

ax1.plot(x,y,linestyle=" ",marker="x")
ax1.plot(Xx,Yx,linestyle=" ",marker="x")
ax1.plot(Xy,Yy,linestyle=" ",marker="x")

# Function for G(x,y)
#-------------------------------------------------- 
nx = 1
ny = 2
n  = nx + ny + 1      # No of free parameters/Conditions to satisfy

m = 2       # No of incident points
x = zeros(Float64,m)
y = zeros(Float64,m)

# First incident point (Fixed Point of the system)
if m > 0
  x = FAC*[0.0  12.0]
  y = FAC*[0.0  4.0]
end  

# X-Derivative Points
mx = 0
Xx = zeros(Float64,mx)
Yx = zeros(Float64,mx)
if (mx>0)
  Xx = x[2]
  Yx = y[2] 
end  

# Y-Derivative Points
my = 1
Xy = zeros(Float64,my)
Yy = zeros(Float64,my)

if my>0
  Xy[1] = x[2]
  Yy[1] = y[2]
end  

gc = NullClineFcn(x,y,Xx,Yx,Xy,Yy,nx,ny)

#c = M\rhs
gc0 = gc[1]
gcx = gc[2:nx+1]
gcy = gc[nx+2:nx+ny+1]

g(x,y) = FXY(x,y,gc0,gcx,gcy)

xi = -10.0
yr0 = -50.0
yr1 = 1.0
dτ  = 1.0e-3
nsteps = 100000

xx,yy = NullClines(g,xi,yr0,yr1,nsteps,dτ)

ax1.plot(xx,yy)

ax1.plot(x,y,linestyle=" ",marker="s")
ax1.plot(Xx,Yx,linestyle=" ",marker="s")
ax1.plot(Xy,Yy,linestyle=" ",marker="s")

grid("on")

pause(0.01)

println("Press x to stop. Any other key to continue")
xin = readline()
if xin !="x"
 
#  close("all")

  ϵ           = 0.1
  if f(0.0,20.0)>0
    F(x,y) = -f(x,y)/ϵ
  else
    F(x,y) =  f(x,y)/ϵ
  end  

  η  = 0.1
  if g(20.0,0.0)>0
    G(x,y) = -η*g(x,y)
  else
    G(x,y) =  η*g(x,y)
  end
  Flow(x,y) = [G(x,y) F(x,y)]
  include("time_stepper_multiple.jl")
end




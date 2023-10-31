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
    s = s + cx[i]*x^i
  end

  for j in 1:ny
    s = s + cy[j]*y^j
  end  

  return s
end  
#---------------------------------------------------------------------- 
function FX(x,c0,cx)
  nx = length(cx)

  s = c0
  for i in 1:nx
    s = s + cx[i]*x^i
  end

  return s
end  
#---------------------------------------------------------------------- 
function FY(y,c0,cy)
  ny = length(cy)

  s = c0
  for i in 1:ny
    s = s + cy[i]*y^i
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

close("all")

lafs = 16

nx = 1
ny = 3
n  = nx + ny + 1      # No of free parameters/Conditions to satisfy

m = 3       # No of incident points
x = zeros(Float64,m)
y = zeros(Float64,m)
if m > 0
  x     = [0.0 -5.0 -5.0]
  y     = [0.0  1.0  6.0]
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

h1 = figure(num=1)
plot(xx,yy)
ax1 = h1.axes[1]
ax1.set_xlabel(L"x", fontsize=lafs)
ax1.set_ylabel(L"y", fontsize=lafs)

plot(x,y,linestyle=" ",marker="x")
plot(Xx,Yx,linestyle=" ",marker="x")
plot(Xy,Yy,linestyle=" ",marker="x")

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
  x = [0.0  8.0]
  y = [0.0  4.0]
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

h1 = figure(num=1)
plot(xx,yy)
ax1 = h1.axes[1]
#ax1.set_xlabel(L"x", fontsize=lafs)
#ax1.set_ylabel(L"y", fontsize=lafs)

plot(x,y,linestyle=" ",marker="s")
plot(Xx,Yx,linestyle=" ",marker="s")
plot(Xy,Yy,linestyle=" ",marker="s")


grid("on")

#ϵ      = 0.1
#H(x,y) = [f(x,y)/ϵ g(x,y)]







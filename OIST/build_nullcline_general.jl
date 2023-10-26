#!/bin/julia

println("Building a null cline based on given Points/gradients by Minimizing the Lagrangian (Error)")


using LinearAlgebra
using Roots
using PyPlot
using Optim

include("NullClines.jl")
include("ElementOf.jl")

function FXY(x,y,c)

  n  = length(x)

  s0 = c[1]*(x.^0)
  s1 = c[2]*x     .+   c[3]*y
  s2 = c[4]*x.^2  .+   c[5]*x.*y       .+   c[6]*y.^2 
  s3 = c[7]*x.^3  .+   c[8]*(x.^2).*y  .+   c[9]*x.*(y.^2)       .+   c[10]*y.^3
  s4 = c[11]*x.^4 .+   c[12]*(x.^3).*y .+   c[13]*(x.^2).*(y.^2) .+   c[14]*x.*(y.^3) .+  c[15]*y.^4

  s  = s0 .+ s1 .+ s2 .+ s3 .+ 0*s4

  return s
end  
#---------------------------------------------------------------------- 
function DxF(x,y,c)

  n  = length(x)

  s0 = 0.0*c[1]*(x.^0)
  s1 = 1*c[2]*x.^0  .+ 0*c[3]*y
  s2 = 2*c[4]*x.^1  .+ 1*c[5]*y          .+ 0*c[6]*y.^2 
  s3 = 3*c[7]*x.^2  .+ 2*c[8]*(x.^1).*y  .+ 1*c[9]*(y.^2)          .+ c[10]*y.^3
  s4 = 4*c[11]*x.^3 .+ 3*c[12]*(x.^2).*y .+ 2*c[13]*(x.^1).*(y.^2) .+ 1*c[14]*(y.^3) .+ 0*c[15]*y.^4

  ds  = s0 .+ s1 .+ s2 .+ s3 .+ 0*s4

  return ds
end  
#---------------------------------------------------------------------- 
function DyF(x,y,c)

  n  = length(y)

  s0 = 0.0*c[1]*(y.^0)
  s1 = 0*c[2] *x.^1  .+ 1*c[3]*(y.^0)
  s2 = 0*c[4] *x.^2  .+ 1*c[5]*x        .+ 2*c[6] *y 
  s3 = 0*c[7] *x.^3  .+ 1*c[8]*(x.^2)   .+ 2*c[9] *x.*y        .+ 3*c[10]*y.^2
  s4 = 0*c[11]*x.^4  .+ 1*c[12]*(x.^3)  .+ 2*c[13]*(x.^2).*y   .+ 3*c[14]*x.*(y.^2) .+ 4*c[15]*y.^3

  ds  = s0 .+ s1 .+ s2 .+ s3 .+ 0*s4

  return ds
end  
#---------------------------------------------------------------------- 
function TERR_F(x,y,xdx,ydx,xdy,ydy,c)

  f   = FXY(x,y,c)
  Er1 = f'*f
  fx  = DxF(xdx,ydx,c)
  Er2 = fx'*fx
  fy  = DyF(xdy,ydy,c)
  Er3 = fy'*fy

  er = Er1 + Er2 + Er3 

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

# First incident point (Fixed Point of the system)
x[1] = 1.0
y[1] = 1.0

# Second incident point
x[2] = -3.0
y[2] = 1.5

# Third incident point 
x[3] = -10.0
y[3] = 5.0

## Fourth incident point 
#x[4] = -10.0
#y[4] = 3.00

## Fifth incident point 
#x[5] = -30.0
#y[5] = 3.5


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


#  Xy[2] = x[3]
#  Yy[2] = y[3]
 
end  

TEC(c) = TERR_F(x,y,Xx,Yx,Xy,Yy,c)

cc = rand(15)

result = optimize(TEC, cc)
c = result.minimizer

f(x,y) = FXY(x,y,c)
#
xi = -10.0
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
#
plot(x,y,linestyle=" ",marker="o")
plot(Xx,Yx,linestyle=" ",marker="o")
plot(Xy,Yy,linestyle=" ",marker="o")
#
#
#fx(x) = FX(x,c0,cx)
#fy(y) = FY(y,0.0,cy)
#
#
#grid("on")









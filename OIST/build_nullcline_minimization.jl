#!/bin/julia

println("Building a null cline based on given Points/gradients by Minimizing the Lagrangian (Error)")


using LinearAlgebra
using Roots
using PyPlot

include("NullClines.jl")
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
#
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


# Build the matrix
N = m + mx + my         # No of constraints to satisfy 
Mat = zeros(Float64,n,n)

# Incident conditions
for i in 1:n
  for j in 1:n 
#   Incident Point contributions      
    for k in 1:m
      if i<=nx+1
        deriv    = (x[k]^(i-1))
      elseif i>nx+1
        deriv    = (y[k]^(i-nx-1))
      end  
     
      if j <= (nx+1)
        Mat[i,j] = Mat[i,j] + (x[k]^(j-1))*deriv
      else
        Mat[i,j] = Mat[i,j] + (y[k]^(j-nx-1))*deriv
      end  
    end     # k
  end       # j
end         # i

## Vanishing Gradient contributions df/dx = 0
Matx = zeros(Float64,n,n)
for i in 2:nx+1
  for j in 2:nx+1 
#   x-gradient contributions 
    for k in 1:mx
      deriv     = (i-1)*(Xx[k]^(i-2))
      Matx[i,j] = Matx[i,j] + (j-1)*(Xx[k]^(j-2))*deriv
    end     # k
  end       # j
end         # i

## Vanishing Gradient contributions df/dy = 0
Maty = zeros(Float64,n,n)
for i in nx+2:n
  for j in nx+2:n 
#   y-gradient contributions 
    for k in 1:my
      l         = i-(nx+1)
      h         = j-(nx+1)
      deriv     = l*(Yy[k]^(l-1))
      Maty[i,j] = Maty[i,j] + h*(Yy[k]^(h-1))*deriv
    end     # k
  end       # j
end         # i

M = Mat + Matx + Maty

F = eigen(M)
j = argmin(F.values)
r0 = minimum(F.values)

println("Minimum abs Eigenvalue: $r0")

c = F.vectors[:,j]

#c = M\rhs
c0 = c[1]
cx = c[2:nx+1]
cy = c[nx+2:nx+ny+1]

f(x,y) = FXY(x,y,c0,cx,cy)

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

plot(x,y,linestyle=" ",marker="o")
plot(Xx,Yx,linestyle=" ",marker="o")
plot(Xy,Yy,linestyle=" ",marker="o")


fx(x) = FX(x,c0,cx)
fy(y) = FY(y,0.0,cy)

#xin = -2.0:0.001:2
#fout = zeros(Float64,length(xin))
#for j in 1:length(xin)
#  fout[j] = fx(xin[j])
#end
#
#plot(xin,fout)

#yin = -15.0:0.01:15
#fout = zeros(Float64,length(yin))
#for j in 1:length(yin)
#  fout[j] = FXY(xi,yin[j],c0,cx,cy)
#end
#
#plot(yin,fout)
#ax1.set_xlabel(L"y", fontsize=lafs)
#ax1.set_ylabel(L"f(y)", fontsize=lafs)


grid("on")









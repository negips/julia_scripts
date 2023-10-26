#!/bin/julia

println("Building a null cline based on given Points/gradients")


using LinearAlgebra
using Roots
using PyPlot

include("NullClines.jl")

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

lafs = 16

nx = 2
ny = 1
n  = nx + ny + 1      # No of free parameters/Conditions to satisfy

m = 4       # No of incident points

x = zeros(Float64,n)
y = zeros(Float64,n)

# First incident point (Fixed Point of the system)
x[1] = 1.0
y[1] = 0.0

# Second incident point
dx = -2.0
dy = -0.0

x[2] = x[1] + dx
y[2] = y[1] + dy

# Third incident point 
x[3] = 0.0
y[3] = 3.0

# Third incident point 
x[4] = 0.5
y[4] = 2.0


# dx/dy = 0.0 at (x3,y3)
x3 = x[3]
y3 = y[3]

# dx/dy = 0.0 at (x4,y4)
x4 = x[3]
y4 = y[3]


# Build the matrix

M = zeros(Float64,n,n)

# Incident conditions
for i in 1:m
  
  # [x[i] y[i] y[i]^2 y[i]^3]
  # x
  M[i,1] = 1.0

  for j in 1:nx
    M[i,j+1] = x[i]^j
  end

  # y
  for j in 1:ny
    k = nx+1+j
    M[i,k] = y[i]^j
  end

end  

## Vanishing Gradients df/dx = 0
#M[m+1,1] = 0.0
#for i in 1:nx
#  M[m+1,i+1] = (i)*x3^(i-1)
#end  

#M[m+1,:] = [0.0 1.0*(y3^0) 2.0*(y3^1) 3.0*(y3^2)]
##M[m+2,:] = [0.0 1.0*(y4^0) 2.0*(y4^1) 3.0*(y4^2)]

rhs = zeros(Float64,n)
for i in 1:m
  rhs[i] = 0.0
end  

AA = qr(M)

r = abs.(diag(AA.R))

j = argmin(r)
r0 = minimum(r)

println("Minimum abs Eigenvalue: $r0")

c = AA.Q[:,j]

#c = M\rhs
c0 = c[1]
cx = c[2:nx+1]
cy = c[nx+2:nx+ny+1]

f(x,y) = FXY(x,y,c0,cx,cy)


# xi = 1000.0
# yr0 = -20.0
# yr1 = 20.0
# dτ  = -1.0e-2
# nsteps = 100000
# 
# xx,yy = NullClines(f,xi,yr0,yr1,nsteps,dτ)
# 
# h1 = figure(num=1)
# plot(xx,yy)
# ax1 = h1.axes[1]
# ax1.set_xlabel(L"x", fontsize=lafs)
# ax1.set_ylabel(L"y", fontsize=lafs)

y = 0.0
fx(x) = f(x,y)

xin = -10.0:0.01:10

fout = zeros(Float64,length(xin))
for j in 1:length(xin)
  fout[j] = fx(xin[j])
end

plot(xin,fout)
grid()









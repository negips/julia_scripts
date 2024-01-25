#!/bin/julia

println("Get Multidimensional rst coordinates of a point")

using PolynomialBases
using LinearAlgebra
using PyPlot

# Performing a Newton Method for approximation

# In 2D
function f2d(r::Float64,s::Float64,x::Vector{Float64},y::Vector{Float64})

  xi = x[1]*(1.0-r)*(1.0-s) + x[2]*r*(1.0-s) + 
       x[3]*(1.0-r)*s       + x[4]*r*s

  yi = y[1]*(1.0-r)*(1.0-s) + y[2]*r*(1.0-s) +
       y[3]*(1.0-r)*s       + y[4]*r*s

 return xi,yi
end  
#---------------------------------------------------------------------- 
# In 3D
function f3d(r::Float64,s::Float64,t::Float64,x::Vector{Float64},y::Vector{Float64},z::Vector{Float64})

  @assert length(x) == length(y) == length(z) "Unequal Lengths of x,y,z"

  xi = x[1]*(1.0-r)*(1.0-s)*(1.0-t) + x[2]*r*(1.0-s)*(1.0-t) + 
       x[3]*(1.0-r)*s*(1.0-t)       + x[4]*r*s*(1.0-t)       +
       x[5]*(1.0-r)*(1.0-s)*t       + x[6]*r*(1.0-s)*t       +
       x[7]*(1.0-r)*s*t             + x[8]*r*s*t

  yi = y[1]*(1.0-r)*(1.0-s)*(1.0-t) + y[2]*r*(1.0-s)*(1.0-t) + 
       y[3]*(1.0-r)*s*(1.0-t)       + y[4]*r*s*(1.0-t)       +
       y[5]*(1.0-r)*(1.0-s)*t       + y[6]*r*(1.0-s)*t       +
       y[7]*(1.0-r)*s*t             + y[8]*r*s*t

  zi = z[1]*(1.0-r)*(1.0-s)*(1.0-t) + z[2]*r*(1.0-s)*(1.0-t) + 
       z[3]*(1.0-r)*s*(1.0-t)       + z[4]*r*s*(1.0-t)       +
       z[5]*(1.0-r)*(1.0-s)*t       + z[6]*r*(1.0-s)*t       +
       z[7]*(1.0-r)*s*t             + z[8]*r*s*t

 return xi,yi,zi
end  
#---------------------------------------------------------------------- 
function Gradf2d(r::Float64,s::Float64,x::Vector{Float64},y::Vector{Float64})

  @assert length(x) == length(y) == length(z) "Unequal Lengths of x,y,z"

  Gf      = zeros(Float64,2,2)
  Gf[1,1] = x[1]*(-1.0)*(1.0-s) + x[2]*(1.0)*(1.0-s) + 
            x[3]*(-1.0)*s       + x[4]*(1.0)*s

  Gf[1,2] = x[1]*(1.0-r)*(-1.0) + x[2]*r*(-1.0) + 
            x[3]*(1.0-r)*(1.0)  + x[4]*r*(1.0)

  Gf[2,1] = y[1]*(-1.0)*(1.0-s) + y[2]*(1.0)*(1.0-s) + 
            y[3]*(-1.0)*s       + y[4]*(1.0)*s

  Gf[2,2] = y[1]*(1.0-r)*(-1.0) + y[2]*r*(-1.0) + 
            y[3]*(1.0-r)*(1.0)  + y[4]*r*(1.0)

 return Gf
end  
#---------------------------------------------------------------------- 
function Gradf3d(r::Float64,s::Float64,t::Float64,x::Vector{Float64},y::Vector{Float64},z::Vector{Float64})

  Gf      = zeros(Float64,3,3)
  Gf[1,1] = x[1]*(-1.0)*(1.0-s)*(1.0-t) + x[2]*(1.0)*(1.0-s)*(1.0-t) + 
            x[3]*(-1.0)*s*(1.0-t)       + x[4]*(1.0)*s*(1.0-t)       +
            x[5]*(-1.0)*(1.0-s)*t       + x[6]*(1.0)*(1.0-s)*t       +
            x[7]*(-1.0)*s*t             + x[8]*(1.0)*s*t

  Gf[1,2] = x[1]*(1.0-r)*(-1.0)*(1.0-t) + x[2]*r*(-1.0)*(1.0-t) + 
            x[3]*(1.0-r)*(1.0)*(1.0-t)  + x[4]*r*(1.0)*(1.0-t)  +
            x[5]*(1.0-r)*(-1.0)*t       + x[6]*r*(-1.0)*t       +
            x[7]*(1.0-r)*(1.0)*t        + x[8]*r*(1.0)*t

  Gf[1,3] = x[1]*(1.0-r)*(1.0-s)*(-1.0) + x[2]*r*(1.0-s)*(-1.0) + 
            x[3]*(1.0-r)*s*(-1.0)       + x[4]*r*s*(-1.0)       +
            x[5]*(1.0-r)*(1.0-s)*(1.0)  + x[6]*r*(1.0-s)*(1.0)  +
            x[7]*(1.0-r)*s*(1.0)        + x[8]*r*s*(1.0)

  Gf[2,1] = y[1]*(-1.0)*(1.0-s)*(1.0-t) + y[2]*(1.0)*(1.0-s)*(1.0-t) + 
            y[3]*(-1.0)*s*(1.0-t)       + y[4]*(1.0)*s*(1.0-t)       +
            y[5]*(-1.0)*(1.0-s)*t       + y[6]*(1.0)*(1.0-s)*t       +
            y[7]*(-1.0)*s*t             + y[8]*(1.0)*s*t

  Gf[2,2] = y[1]*(1.0-r)*(-1.0)*(1.0-t) + y[2]*r*(-1.0)*(1.0-t) + 
            y[3]*(1.0-r)*(1.0)*(1.0-t)  + y[4]*r*(1.0)*(1.0-t)  +
            y[5]*(1.0-r)*(-1.0)*t       + y[6]*r*(-1.0)*t       +
            y[7]*(1.0-r)*(1.0)*t        + y[8]*r*(1.0)*t

  Gf[2,3] = y[1]*(1.0-r)*(1.0-s)*(-1.0) + y[2]*r*(1.0-s)*(-1.0) + 
            y[3]*(1.0-r)*s*(-1.0)       + y[4]*r*s*(-1.0)       +
            y[5]*(1.0-r)*(1.0-s)*(1.0)  + y[6]*r*(1.0-s)*(1.0)  +
            y[7]*(1.0-r)*s*(1.0)        + y[8]*r*s*(1.0)


  Gf[3,1] = z[1]*(-1.0)*(1.0-s)*(1.0-t) + z[2]*(1.0)*(1.0-s)*(1.0-t) + 
            z[3]*(-1.0)*s*(1.0-t)       + z[4]*(1.0)*s*(1.0-t)       +
            z[5]*(-1.0)*(1.0-s)*t       + z[6]*(1.0)*(1.0-s)*t       +
            z[7]*(-1.0)*s*t             + z[8]*(1.0)*s*t

  Gf[3,2] = z[1]*(1.0-r)*(-1.0)*(1.0-t) + z[2]*r*(-1.0)*(1.0-t) + 
            z[3]*(1.0-r)*(1.0)*(1.0-t)  + z[4]*r*(1.0)*(1.0-t)  +
            z[5]*(1.0-r)*(-1.0)*t       + z[6]*r*(-1.0)*t       +
            z[7]*(1.0-r)*(1.0)*t        + z[8]*r*(1.0)*t

  Gf[3,3] = z[1]*(1.0-r)*(1.0-s)*(-1.0) + z[2]*r*(1.0-s)*(-1.0) + 
            z[3]*(1.0-r)*s*(-1.0)       + z[4]*r*s*(-1.0)       +
            z[5]*(1.0-r)*(1.0-s)*(1.0)  + z[6]*r*(1.0-s)*(1.0)  +
            z[7]*(1.0-r)*s*(1.0)        + z[8]*r*s*(1.0)
           
 return Gf
end  
#---------------------------------------------------------------------- 

function Get_rs_Newton(x0::Float64,y0::Float64,x::Vector{Float64},y::Vector{Float64})

  tol   = 5.0*eps()
  maxit = 20

  # Initial guesses
  r = 0.5 
  s = 0.5

  xi,yi = f2d(r,s,x,y)

  dx = x0-xi
  dy = y0-yi

  if abs(dx)>0.0 || abs(dy)>0.0 
    if abs(dy)>(dx)
      res = abs(dy)*sqrt((dx/dy)^2 + 1.0)
    else
      res = abs(dx)*sqrt((dy/dx)^2 + 1.0)
    end
  else
    res = 0.0
  end  
    
  it = 0
  while it<maxit && res>tol
    it     += 1
    GF      = Gradf2d(r,s,x,y)

    drs = -inv(GF)*[xi-x0;
                    yi-y0]
    r   = r + drs[1]
    s   = s + drs[2]

    xi,yi   = f2d(r,s,x,y)

    dx = x0-xi
    dy = y0-yi

    if abs(dy)>0.0 || abs(dx)>0.0
      if abs(dy)>abs(dx)
        res = abs(dy)*sqrt((dx/dy)^2 + 1.0)
      else
        res = abs(dx)*sqrt((dy/dx)^2 + 1.0)
      end
    else
      res = 0.0
    end  
  end       # while

  #println("Iterations: ", it)
  #println("Residual: ", res)

  return r,s
end  
#----------------------------------------------------------------------
function Get_rst_Newton(x0::Float64,y0::Float64,z0::Float64,x::Vector{Float64},y::Vector{Float64},z::Vector{Float64})

  tol   = 5.0*eps()
  maxit = 30

  # Initial guesses
  r = 0.5 
  s = 0.5
  t = 0.5

  xi,yi,zi = f3d(r,s,t,x,y,z)

  dx = x0-xi
  dy = y0-yi
  dz = z0-zi

  dm = max(abs(dx),abs(dy),abs(dz))
  if dm>0.0
    res = dm*sqrt((dx/dm)^2 + (dy/dm)^2 + (dz/dm)^2)
  else
    res = 0.0
  end
    
  it = 0
  while it<maxit && res>tol
    it     += 1
    GF      = Gradf3d(r,s,t,x,y,z)

    drst    = -inv(GF)*[xi-x0;
                        yi-y0;
                        zi-z0]

    r   = r + drst[1]
    s   = s + drst[2]
    t   = t + drst[3]

    xi,yi,zi   = f3d(r,s,t,x,y,z)

    dx = x0-xi
    dy = y0-yi
    dz = z0-zi

    dm = max(abs(dx),abs(dy),abs(dz))
    if dm>0.0
      res = dm*sqrt((dx/dm)^2 + (dy/dm)^2 + (dz/dm)^2)
    else
      res = 0.0
    end
   
  end       # while

  #println("Iterations: ", it)
  #println("Residual: ", res)

  return r,s,t
end  
#----------------------------------------------------------------------

close("all")

dim  = 2
n    = 2^dim

x2   = zeros(Float64,n)
y2   = zeros(Float64,n)

x2[1] = 0.0
x2[2] = 4.0
x2[3] = 2.5
x2[4] = 5.0

y2[1] = 0.0
y2[2] = 0.0
y2[3] = 2.0
y2[4] = 4.0

plot(x2,y2, linestyle="none",marker="o") 

r0   = 0.01
s0   = 0.94

xi, yi = f2d(r0,s0,x2,y2)
plot(xi,yi,linestyle="none",marker="s")

r,s    = Get_rs_Newton(xi,yi,x2,y2)

println("r,s: ",r,",",s)

# 3D test
dim  = 3
n    = 2^dim

x3 = zeros(Float64,n)
y3 = zeros(Float64,n)
z3 = zeros(Float64,n)

x3[1] = 0.0
x3[2] = 4.0
x3[3] = 2.5
x3[4] = 5.0
x3[5] = 1.0
x3[6] = 5.0
x3[7] = 3.5
x3[8] = 6.0


y3[1] = 0.0
y3[2] = 0.0
y3[3] = 2.0
y3[4] = 4.0
y3[5] = 1.0
y3[6] = 1.0
y3[7] = 3.0
y3[8] = 5.0


z3[1] = 0.0
z3[2] = 0.0
z3[3] = 0.0
z3[4] = 0.0
z3[5] = 1.0
z3[6] = 1.0
z3[7] = 2.5
z3[8] = 2.5


r0   = 0.91
s0   = 0.94
t0   = 0.01

xi, yi, zi = f3d(r0,s0,t0,x3,y3,z3)

h2  = figure(num=2)
plot3D(x3,y3,z3,linestyle="none",marker="o")
ax = h2.gca()
ax.set_xlabel("X")
ax.set_ylabel("Y")
ax.set_zlabel("Z")

r,s,t    = Get_rst_Newton(xi,yi,zi,x3,y3,z3)

println("r,s,t: ",r,",",s,",",t)







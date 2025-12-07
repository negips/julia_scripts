#!/bin/julia

println("Langmuir-Hinshelwood mechanism in chemical kinetics.")

using LinearAlgebra
using PyPlot


#---------------------------------------------------------------------- 
function timederivative(x::Float64,y::Float64,k1::Float64,k2::Float64,k_2::Float64,k3::Float64)

  xdot = 2.0*k1*(1.0 - x - y)^2         - k3*x*y
  ydot = 1.0*k2*(1.0 - x - y)   - k_2*y - k3*x*y

  return xdot,ydot
end 
#---------------------------------------------------------------------- 
# 4th Order Runge-Kutta Steps
function RK4(x::Float64,y::Float64,dt::Float64,TD)

  two = 2.0
  six = 6.0

  xd1,yd1 = TD(x,y) 
  x1 = x  + dt/two*xd1
  y1 = y  + dt/two*yd1

  xd2,yd2 = TD(x1,y1) 
  x2 = x  + dt/two*xd2
  y2 = y  + dt/two*yd2

  xd3,yd3 = TD(x2,y2) 
  x3 = x  + dt*xd3
  y3 = y  + dt*yd3

  xd4,yd4 = TD(x3,y3) 
  xo = x  + dt/six*(xd1 + two*xd2 + two*xd3 + xd4)
  yo = y  + dt/six*(yd1 + two*yd2 + two*yd3 + yd4)

  return xo,yo
end  
#---------------------------------------------------------------------- 


k3    = 10.0
k1    = 0.1
k2    = 0.1
k_2   = 0.1
dt    = 0.01

TD(x,y)    = timederivative(x,y,k1,k2,k_2,k3)

x     = 1.0
y     = 1.1

nsteps = 1000000
xhist  = zeros(nsteps)
yhist  = zeros(nsteps)

for i in 1:nsteps
  global x,y

  x,y = RK4(x,y,dt,TD)

  xhist[i] = x
  yhist[i] = y

end  

# close("all")

 plot(xhist,yhist)









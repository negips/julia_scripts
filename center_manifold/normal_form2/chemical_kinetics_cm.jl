#!/bin/julia

println("Langmuir-Hinshelwood mechanism in chemical kinetics.")

using LinearAlgebra
using PyPlot


#---------------------------------------------------------------------- 
function deviation_TD(x::Float64,y::Float64,μ1::Float64,μ2::Float64,μ3::Float64)

  local x0   = 1.0
  local y0   = 0.0

  xdot = 2.0*μ1*(1.0 - (x0 + x) - (y0 + y))^2 - (x0 + x)*(y0 + y)
  ydot = 1.0*μ2*(1.0 - (x0 + x) - (y0 + y))   - μ3*(y0 + y)  - (x0 + x)*(y0 + y)

  return xdot,ydot
end 
#---------------------------------------------------------------------- 
function deviation_TD2(x::Float64,y::Float64,μ1::Float64,μ2::Float64,μ3::Float64)

  xdot =  2.0*μ1*(x + y)^2         - (1.0*y + x*y)
  ydot = -1.0*μ2*(x + y)   - μ3*y  - (1.0*y + x*y)

  return xdot,ydot
end 
#---------------------------------------------------------------------- 

# 4th Order Runge-Kutta Steps
function RK4(xi::Float64,yi::Float64,dt::Float64,TD)

  two = 2.0
  six = 6.0

  xd0,yd0 = TD(xi,yi) 
  x1 = xi  + dt/two*xd0
  y1 = yi  + dt/two*yd0

  xd1,yd1 = TD(x1,y1) 
  x2 = xi  + dt/two*xd1
  y2 = yi  + dt/two*yd1

  xd2,yd2 = TD(x2,y2) 
  x3 = xi  + dt*xd2
  y3 = yi  + dt*yd2

  xd3,yd3 = TD(x3,y3) 
  xo = xi  + dt/six*(xd0 + two*xd1 + two*xd2 + xd3)
  yo = yi  + dt/six*(yd0 + two*yd1 + two*yd2 + yd3)

  return xo,yo
end  
#---------------------------------------------------------------------- 


μ1    = 0.1
μ2    = 0.1
μ3    = 0.1
dt    = 0.01

# Deviation Time Derivative
TD(x,y)    = deviation_TD(x,y,μ1,μ2,μ3)
TD2(x,y)   = deviation_TD2(x,y,μ1,μ2,μ3)

xi     = 0.1
yi     = 0.1

nsteps = 4000000
xhist  = zeros(nsteps)
yhist  = zeros(nsteps)

for i in 1:nsteps
  global xi,yi

  xi,yi = RK4(xi,yi,dt,TD2)

  xhist[i] = xi
  yhist[i] = yi

end  

close("all")

plot(xhist,yhist)









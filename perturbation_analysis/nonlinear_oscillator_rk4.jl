println("Nonlinear Oscillator")
println("Model: x'' + αx' + x + ϵx^3 = 0")

import Pkg
#Pkg.add("PyPlot")
#Pkg.add("Blink")

using Blink
using PyPlot,Colors,PyCall
using LinearAlgebra
using Random

#---------------------------------------------------------------------- 
function nonosci_rk4!(y,x,α,ϵ,dt)

  dy  = -(α*y + x + ϵ*x^3)
  dx  = y
  y1  = y + dt/2.0*dy
  x1  = x + dt/2.0*dx

  dy1 = -(α*y1 + x1 + ϵ*x1^3)
  dx1 = y1
  y2  = y + dt/2.0*dy1
  x2  = x + dt/2.0*dx1

  dy2 = -(α*y2 + x2 + ϵ*x2^3)
  dx2 = y2
  y3  = y + dt*dy2
  x3  = x + dt*dx2

  dy3 = -(α*y3 + x3 + ϵ*x3^3)
  dx3 = y3
  y   = y + dt/6.0*(dy + 2.0*dy1 + 2.0*dy2 + dy3)
  x   = x + dt/6.0*(dx + 2.0*dx1 + 2.0*dx2 + dx3)

  return y,x
end
#---------------------------------------------------------------------- 


# Initial condition
y0 = 2.0                # dx/dt
x0 = 0.0                # x

y  = copy(y0)
x  = copy(x0)

α = -0.0          # Negative is unstable
ϵ = 0.2

nsteps = 50000
dt = 0.01

xhist    = zeros(Float64,nsteps+1)
yhist    = zeros(Float64,nsteps+1)

xhist[1] = x0 
yhist[1] = y0 

for i in 1:nsteps
  global x,y
  global xhist,yhist

  y,x = nonosci_rk4!(y,x,α,ϵ,dt)

  yhist[i+1] = y
  xhist[i+1] = x
end

time = range(0.0,step=dt,length=nsteps+1)./π
plot(time,xhist);
ax = gca();
ax.grid("on")

#---------------------------------------------------------------------- 






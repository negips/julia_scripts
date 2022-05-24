println("Rayleigh Oscillator")
println("Model: x'' + αx' + (ω^2)*x - ϵ*x' + ϵ/3.0*(x')^3 = 0")

import Pkg
#Pkg.add("PyPlot")
#Pkg.add("Blink")

using Blink
using PyPlot,Colors,PyCall
using LinearAlgebra
using Random

#---------------------------------------------------------------------- 
function rayosci_rk4!(y,x,α,ω,ϵ,dt)

  dy  = -(α*y + (ω^2)*x + -ϵ*y + ϵ/3.0*y^3)
  dx  = y
  y1  = y + dt/2.0*dy
  x1  = x + dt/2.0*dx

  dy1 = -(α*y1 + (ω^2)*x1 -ϵ*y1 + ϵ/3.0*y1^3)
  dx1 = y1
  y2  = y + dt/2.0*dy1
  x2  = x + dt/2.0*dx1

  dy2 = -(α*y2 + (ω^2)*x2 -ϵ*y2 + ϵ/3.0*y2^3)
  dx2 = y2
  y3  = y + dt*dy2
  x3  = x + dt*dx2

  dy3 = -(α*y3 + (ω^2)*x3 -ϵ*y3 + ϵ/3.0*y3^3)
  dx3 = y3
  y   = y + dt/6.0*(dy + 2.0*dy1 + 2.0*dy2 + dy3)
  x   = x + dt/6.0*(dx + 2.0*dx1 + 2.0*dx2 + dx3)

  return y,x
end
#---------------------------------------------------------------------- 

function modosci_rk4!(y,x,α,ω,ϵ,dt)

  dy  = -(α*y + (ω^2)*x + -0.0*ϵ*y + ϵ/3.0*y^3)
  dx  = y
  y1  = y + dt/2.0*dy
  x1  = x + dt/2.0*dx

  dy1 = -(α*y1 + (ω^2)*x1 -0.0*ϵ*y1 + ϵ/3.0*y1^3)
  dx1 = y1
  y2  = y + dt/2.0*dy1
  x2  = x + dt/2.0*dx1

  dy2 = -(α*y2 + (ω^2)*x2 -0.0*ϵ*y2 + ϵ/3.0*y2^3)
  dx2 = y2
  y3  = y + dt*dy2
  x3  = x + dt*dx2

  dy3 = -(α*y3 + (ω^2)*x3 -0.0*ϵ*y3 + ϵ/3.0*y3^3)
  dx3 = y3
  y   = y + dt/6.0*(dy + 2.0*dy1 + 2.0*dy2 + dy3)
  x   = x + dt/6.0*(dx + 2.0*dx1 + 2.0*dx2 + dx3)

  return y,x
end
#---------------------------------------------------------------------- 

# Initial condition
a  = 0.05
y0 = 2.0*a
x0 = 0.0

y  = y0
x  = x0

α = -0.00          # Negative is unstable
ω = 1.2
ϵ = 0.2

nsteps = 1000000
dt = 0.0001

yhist = zeros(Float64,nsteps+1)
xhist = zeros(Float64,nsteps+1)

yhist[1] = y0
xhist[1] = x0

for i in 1:nsteps
  global x,y
  global xhist,yhist

  y,x = rayosci_rk4!(y,x,α,ω,ϵ,dt)
#  y,x = modosci_rk4!(y,x,α,ω,ϵ,dt)

  xhist[i+1] = x
  yhist[i+1] = y
end

t = range(0.0,step=dt,length=nsteps+1)
plot(t,xhist);

# first order perturbation solution

Ω  = ω*sqrt(Complex(1 - (σ/ω)^2))
σ  = -α/2.0
λ  = σ + im*Ω
β  = λ*exp.(σ*t)
β2 = (λ*λ').*exp.(2.0*σ*t)
θ  = -π/2.0
R0 = a/(abs(λ))
γi = 0.5*α/Ω

A0     = exp.(λ*t)
R      = 2.0*R0*(exp.(-ϵ*t) .+ (R0^2).*(β2).*(1.0 .- exp.(-ϵ*t))).^(-1/2)
θ      = -γi*t .+ γi*log.(β2*R0^2.0.*(exp.(t) .- 1.0) .+ 1.0) .- π/2 
expθ   = exp.(im*θ)
zpert1 = R.*expθ.*A0
xpert1 = real.(zpert1)
plot(t,xpert1);






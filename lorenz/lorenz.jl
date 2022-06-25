#!/usr/bin/julia

println("Visualizing the Lorenz attractor")


# ̇x = σ(y - x)
# ̇y = ρx -y -zx
# ̇z = -βz + xy

# Standard parameters:
# σ = 10.0
# ρ = 28.0        # also R
# β = 8.0/3.0

using LinearAlgebra
using PyPlot,Blink

include("Lorenz_ddt.jl")
include("RK4_lorenz.jl")

close("all")

lafs = 16

σ = 10.0
ρ = 28.0          # also R
β = 8.0/3.0

dt = 0.001

nsteps = 10000

x = zeros(Float64,nsteps)
y = zeros(Float64,nsteps)
z = zeros(Float64,nsteps)

xi = 10.0
yi = 10.0
zi = 10.0

for i in 1:nsteps
   
   global xi,yi,zi
   
   xi,yi,zi = RK4_lorenz!(xi,yi,zi,σ,ρ,β,dt)

   x[i] = xi
   y[i] = yi
   z[i] = zi

end   

plot(x,z)
ax1 = gca();
ax1.set_xlabel(L"x",fontsize=lafs)
ax1.set_ylabel(L"z",fontsize=lafs)

println("Done.")







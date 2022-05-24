println("Comparing perturbed oscillation secular series")
println("Model: x'' + (1+ϵ)^2*x = 0")


using PyPlot,Colors,PyCall
using LinearAlgebra

ϵ = 0.01
A = 1.0

dt = 0.01
nsteps = 1000
t = range(0.0,step=dt,length=nsteps)
tc = complex.(t)
ω = 1+ϵ
z = A*exp.(im*ω*t)
x = real.(z)

ω0 = 1.0
z0 = A*exp.(im*ω0*t)
z1 = -2.0*ϵ*t
z2 = ϵ^2*(2.0*t.^2 .- t)
z3 = ϵ^3*(2.0*t.^2 .- 4.0/3.0*t.^3)

zp0 = z0
zp1 = z0.*(1.0 .+ z1)
zp2 = z0.*(1.0 .+ z1 .+ z2)
zp3 = z0.*(1.0 .+ z1 .+ z2 .+ z3)

xp0 = real.(zp0) 
xp1 = real.(zp1) 
xp2 = real.(zp2) 
xp3 = real.(zp3) 

plot(t,x)
plot(t,xp0)
plot(t,xp1)
plot(t,xp2)
plot(t,xp3)


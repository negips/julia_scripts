#!/bin/julia

println("Calculate the phase flow of the Fitzhugh-Nagumo model")

using LinearAlgebra
using Random
using Printf
using PyPlot,PyCall
using Roots

#include("Meinhardt.jl")
#include("Meinhardt_Nullclines.jl")
include("../NullClines.jl")
include("FitzhughNagumo.jl")
include("$JULIACOMMON/GetBDF.jl")
include("$JULIACOMMON/GetEXT.jl")
include("$JULIACOMMON/MoveFigure.jl")

#-------------------------------------------------- 

function Translate(f,g,a,b,dt)

  anew = a + dt*f(a,b) 
  bnew = b + dt*g(a,b)

  return anew,bnew
end  

#-------------------------------------------------- 


close("all")

lafs = 16

a  = -0.5
b  =  0.5
ϵ  =  0.1


f(x,y)   = FitzhughNagumo(x,y,a,b)[1]./ϵ
g(x,y)   = FitzhughNagumo(x,y,a,b)[2]

nstep = 100000
dt    = 0.001

ntraj = 50

cm    = get_cmap("tab20");

u0    = 6.0*(rand(ntraj) .- 0.5) .- 1.32
v0    = 6.0*(rand(ntraj) .- 0.5) .- 1.66

h1 = figure(num=1)

for n in 1:ntraj
  global u,v

  uhist = zeros(nstep,1)
  vhist = zeros(nstep,1)
 
  jj = mod(n-1,20)
  cmap = cm(jj)
  u = u0[n]
  v = v0[n]
  plot(u,v,linestyle="-",marker="o",color=cmap)
  uhist[1] = u
  vhist[1] = v

  for i in 2:nstep
    u,v      = Translate(f,g,u,v,dt)
    uhist[i] = u
    vhist[i] = v
  end
  
  plot(uhist,vhist,linestyle="--",color=cmap)

end  

ax1 = h1.axes[1]
ax1.set_xlabel(L"u", fontsize=lafs)
ax1.set_ylabel(L"v", fontsize=lafs)
#ax1.legend(fontsize=12)

Uu,Uv,Vu,Vv = FitzhughNagumoNullClines(a,b)


cm2   = get_cmap("tab10");
plot(Uu,Uv,linestyle="-",color=cm2(0),linewidth=2)
plot(Vu,Vv,linestyle="-",color=cm2(1),linewidth=2)

title = "a=$a; b=$b"
ax1.set_title(title, fontsize=lafs)
ax1.set_xlim(-4.0,4.0)
ax1.set_ylim(-4.0,4.0)
MoveFigure(h1,1250,500)
legend()










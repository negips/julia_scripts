#!/bin/julia

println("Calculate the phase flow of the Barkley's Model for pipe flow")

using LinearAlgebra
using Random
using Printf
using PyPlot,PyCall
using Roots

#include("Meinhardt.jl")
#include("Meinhardt_Nullclines.jl")
include("../NullClines.jl")
include("Barkley.jl")
include("../GetBDF.jl")
include("../GetEXT.jl")

function Translate(f,g,a,b,dt)

  anew = a + dt*f(a,b) 
  bnew = b + dt*g(a,b)

  return anew,bnew
end  

#-------------------------------------------------- 


close("all")

lafs = 16

σ  = 0.0
r  = 1.2


g(x,y)   = BarkleyPipe(x,y,r)[1]
f(x,y)   = BarkleyPipe(x,y,r)[2]

nstep = 1000000
dt    = 0.001

ntraj = 3

cm    = get_cmap("tab20");

q0    = -1.49*rand(ntraj) .+ 1.5
u0    = -2.0*rand(ntraj) .+ 2.5

h1 = figure(num=1)

for n in 1:ntraj
  global q,u

  qhist = zeros(nstep,1)
  uhist = zeros(nstep,1)
 
  jj = mod(n-1,20)
  cmap = cm(jj)
  q = q0[n]
  u = u0[n]
  plot(u,q,linestyle="-",marker="o",color=cmap)
  qhist[1] = q
  uhist[1] = u

  for i in 2:nstep
    u,q   = Translate(g,f,u,q,dt)
    qhist[i] = q
    uhist[i] = u
  end
  
  plot(uhist,qhist,linestyle="--",color=cmap)

end  

ax1 = h1.axes[1]
ax1.set_xlabel(L"u", fontsize=lafs)
ax1.set_ylabel(L"q", fontsize=lafs)
#ax1.legend(fontsize=12)

Qq,Qu,Uq,Uu = BarkleyNullClines(r)

plot(Qu,Qq,linestyle="-",color="blue",linewidth=2)
plot(Uu,Uq,linestyle="-",color="red",linewidth=2)

#legend()










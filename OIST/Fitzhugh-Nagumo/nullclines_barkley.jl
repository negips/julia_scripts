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
R  = [0.60]


g(x,y)   = BarkleyPipe(x,y,r,σ)[1]
f(x,y)   = BarkleyPipe(x,y,r,σ)[2]

nphase = length(R)

cm    = get_cmap("tab20");

h1 = figure(num=1)

for n in 1:nphase
 
  j1 = mod(2*(n-1),20)
  cmapQ = cm(j1)
  j2 = mod(2*(n-1)+1,20)
  cmapU = cm(j2)

  r = R[n]

  Qq,Qu,Uq,Uu = BarkleyNullClines(r)
  
  plot(Qu,Qq,linestyle="-",color="blue",linewidth=2)
  plot(Uu,Uq,linestyle="-",color="red",linewidth=2)

end  

ax1 = h1.axes[1]
ax1.set_xlabel(L"u", fontsize=lafs)
ax1.set_ylabel(L"q", fontsize=lafs)

if length(R)==1
  title = "R=$(R[1])"
  ax1.set_title(title, fontsize=lafs)
end  

#ax1.legend(fontsize=12)


#legend()


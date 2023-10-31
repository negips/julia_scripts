#!/bin/julia

println("Calculate the phase flow in H. Meinhardt (1995) The algorithmic beauty of sea shells")

using LinearAlgebra
using Random
using Printf
using PyPlot,PyCall
using Roots

include("Meinhardt.jl")
include("Meinhardt_Nullclines.jl")
include("GetBDF.jl")
include("GetEXT.jl")

function Translate(f,g,a,b,dt)

  anew = a + dt*f(a,b) 
  bnew = b + dt*g(a,b)

  return anew,bnew
end  

#-------------------------------------------------- 


#close("all")

include("test_nullcline_fcn.jl")

lafs = 16

ϵ           = 0.1
if f(0.0,20.0)>0
  Adotxy(x,y) = -f(x,y)/ϵ
else
  Adotxy(x,y) =  f(x,y)/ϵ
end  

if g(20.0,0.0)>0
  Bdotxy(x,y) = -g(x,y)
else
  Bdotxy(x,y) =  g(x,y)
end  


nstep = 500000
dt    = 0.0001

ntraj = 20

cm    = get_cmap("tab20");

a0    =  6.0*(rand(ntraj) .- 0.5) .+ 3.0
b0    =  6.0*(rand(ntraj) .- 0.5) .+ 0.

h1 = figure(num=1)

for n in 1:ntraj
  global a,b

  ahist = zeros(nstep,1)
  bhist = zeros(nstep,1)
 
  jj = mod(n-1,20)
  cmap = cm(jj)
  a = a0[n]
  b = b0[n]
  plot(b,a,linestyle="-",marker="o",color=cmap)
  ahist[1] = a
  bhist[1] = b

  for i in 2:nstep
    b,a   = Translate(Bdotxy,Adotxy,b,a,dt)
    ahist[i] = a
    bhist[i] = b
  end
  
  plot(bhist,ahist,linestyle="--",color=cmap)

end  

ax1 = h1.axes[1]
ax1.set_xlabel(L"b", fontsize=lafs)
ax1.set_ylabel(L"a", fontsize=lafs)
#ax1.legend(fontsize=12)











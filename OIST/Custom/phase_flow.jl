#!/bin/julia

println("Calculate the phase flow in H. Meinhardt (1995) The algorithmic beauty of sea shells")

using LinearAlgebra
using Random
using Printf
using PyPlot,PyCall
using Roots

#include("Meinhardt.jl")
#include("Meinhardt_Nullclines.jl")
include("$JULIACOMMON/GetBDF.jl")
include("$JULIACOMMON/GetEXT.jl")

function Translate(f,g,a,b,dt)

  anew = a + dt*f(a,b) 
  bnew = b + dt*g(a,b)

  return anew,bnew
end  

#-------------------------------------------------- 

nstep = 50000
dt    = 0.001

ntraj = 30

cm    = get_cmap("tab20");

#a0    =  6.0*(rand(ntraj) .- 0.5) .+ 3.0
#b0    =  6.0*(rand(ntraj) .- 0.5) .+ 0.

slope  = 0.5
y_x0   = 1.2
b1     = LinRange(-0.7,4.5,ntraj)
a1     = y_x0 .+ slope*(b1) 

ntraj  = 30
slope  = 0.47
y_x0   = 1.0
b2     = LinRange(-0.5,4.8,ntraj)
a2     = y_x0 .+ slope*(b2) 

ntraj  = 30
slope  = 0.001
y_x0   = 7.5
b3     = LinRange(-2.0,7.5,ntraj)
a3     = y_x0 .+ slope*(b3) 

ntraj  = 60
slope  = 0.001
y_x0   = -2.0
b4     = LinRange(-2.0,7.5,ntraj)
a4     = y_x0 .+ slope*(b4) 


a0     = [a1; a2; a3; a4]
b0     = [b1; b2; b3; b4]

ntr    = length(a0)

for n in 1:ntr
  global a,b

  ahist = zeros(nstep,1)
  bhist = zeros(nstep,1)
 
  jj = mod(n-1,20)
  cmap = cm(jj)
  a = a0[n]
  b = b0[n]
  plot(b,a,linestyle="-",marker=".",color="black")
  ahist[1] = a
  bhist[1] = b

  for i in 2:nstep
    b,a   = Translate(G,F,b,a,dt)
    ahist[i] = a
    bhist[i] = b
  end
  
  plot(bhist,ahist,linestyle=":",color="gray")

end  











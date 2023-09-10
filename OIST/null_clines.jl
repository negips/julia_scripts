#!/bin/julia

println("Calculate the null clines for Eqns 2.1a, 2.1b in H. Meinhardt (1995) The algorithmic beauty of sea shells")

using LinearAlgebra
using Random
using Printf
using PyPlot,PyCall

#-------------------------------------------------- 
function Adot(a,b,s,ra,ba)

  adot = s*(a.^2 ./b .- ba) .- ra*a

  return adot
end

function Bdot(a,b,s,rb,bb)

  bdot = s*a.^2 .- rb*b .+ bb

  return bdot
end

function Translate(a,b,dt,s,ra,ba,rb,bb)

  anew = a + dt*Adot(a,b,s,ra,ba)    
  bnew = b + dt*Bdot(a,b,s,rb,bb)    

  return anew,bnew
end  


#-------------------------------------------------- 

close("all")

ra = 1.0
s  = 1.0
ba = 2.0

ca = ra/(2*s)

rb = 5.0
bb = 5.0
cb = rb/s

x1 = 0:0.01:5

A1  = ca*x1 .+ sqrt.(ca*ca*x1.*x1 .+ ba*x1)
A2  = ca*x1 .- sqrt.(ca*ca*x1.*x1 .+ ba*x1)



y1 = -3.0:0.01:3.0

B   = s/rb*y1.*y1 .+ bb/rb 


h1 = figure(num=1)
plot(x1,A1)
plot(x1,A2)

plot(B,y1)

lafs = 16

ax1 = h1.axes[1]
ax1.set_xlabel(L"b", fontsize=lafs)
ax1.set_ylabel(L"a", fontsize=lafs)
#ax1.legend(fontsize=12)

ai = [4.0 0.0 -3.0]
bi = [2.0 2.0 2.0]
adot1 = Adot(ai,bi,s,ra,ba)
bdot1 = Bdot(ai,bi,s,rb,bb)


nstep = 400
dt    = 0.01

a     = -4.12
b     = 2.00

ahist = zeros(nstep,1)
bhist = zeros(nstep,1)

for i in 1:nstep
  global a,b
  a,b   = Translate(a,b,dt,s,ra,ba,rb,bb)
  ahist[i] = a
  bhist[i] = b
end

ax1.plot(bhist,ahist,linestyle="--")






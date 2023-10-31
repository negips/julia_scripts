#!/bin/julia

println("Calculate the phase flow in H. Meinhardt (1995) The algorithmic beauty of sea shells")

using LinearAlgebra
using Random
using Printf
using PyPlot,PyCall
using Roots

include("Meinhardt.jl")
include("Meinhardt_Nullclines.jl")
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

c    = 0.10
d    = 1.00

EQN = "1987_1"
if (EQN == "2.1")
  F(x,y,s)   = Meinhardt_21(x,y,s)
elseif (EQN == "2.4")
  F(x,y,s)   = Meinhardt_24(x,y,s)
elseif (EQN == "2.5")
  F(x,y,s)   = Meinhardt_25(x,y,s)
elseif (EQN == "3.1")
  F(x,y,s)   = Meinhardt_31(x,y,s)
elseif (EQN == "6.1")
  F(x,y,s)   = Meinhardt_61(x,y,c,d,s)
elseif (EQN == "1987_1")
  F(x,y,s)   = Meinhardt_1987_1(x,y)
elseif (EQN == "1987_2")
  F(x,y,s)   = Meinhardt_1987_2(x,y)
elseif (EQN == "1987_2_branching")
  F(x,y,s)   = Meinhardt_1987_2_branching(x,y,s)
else
  display("$EQN not defined.")
end


s  = 1.00
Adotxy(x,y) = F(x,y,s)[1]
Adotx(x) = Adotxy(x,b)
Adoty(y) = Adotxy(a,y)

Bdotxy(x,y) = F(x,y,s)[2]
Bdotx(x) = Bdotxy(x,b)
Bdoty(y) = Bdotxy(a,y)


nstep = 500000
dt    = 0.001

ntraj = 0

cm    = get_cmap("tab20");

a0    =  5.3*(rand(ntraj) .- 0.5) .+ 3.35
b0    =  0.3*(rand(ntraj) .- 0.5) .+ 0.9

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
    a,b   = Translate(Adotxy,Bdotxy,a,b,dt)
    ahist[i] = a
    bhist[i] = b
  end
  
  plot(bhist,ahist,linestyle="--",color=cmap)

end  

if (ntraj > 0)
  ax1 = h1.axes[1]
  ax1.set_xlabel(L"b", fontsize=lafs)
  ax1.set_ylabel(L"a", fontsize=lafs)
end  
#ax1.legend(fontsize=12)

G(x,y) = F(x,y,s)
Aa,Ab,Ba,Bb = Meinhardt_Nullclines(EQN,G)

n1 = size(Aa,2)
if n1>1
  lab1 = Vector{String}(undef,n1)
  lab1[1] = "∂A/∂t=0"
  for i in 2:n1
    lab1[i] = ""
  end
else
  lab1 = "∂A/∂t=0"
end  

n2 = size(Ba,2)
if n2>1
  lab2 = Vector{String}(undef,n2)
  lab2[1] = "∂B/∂t=0"
  for i in 2:n2
    lab2[i] = ""
  end
else
  lab2 = "∂B/∂t=0"
end  

plot(Ab,Aa,linestyle="-",color="blue",linewidth=2, label=lab1)
plot(Bb,Ba,linestyle="-",color="red",linewidth=2,  label=lab2)

if EQN == "3.1"
  Aa2,Ab2 = find_nullcline(Adotxy,-1.0,-1.01,-0.3,0.0,120000,-1.0e-5);
  plot(Ab2,Aa2,linestyle="-",color="black",linewidth=2)
  Aa3,Ab3 = find_nullcline(Adotxy,-1.0,-1.007,0.0,0.1,30000,-1.0e-4);
  plot(Ab3,Aa3,linestyle="-",color="black",linewidth=2)
elseif EQN == "1987_2_branching"
  Aa2,Ab2 = find_nullcline(Adotxy,-2.0,-2.0,-5.0,1.0,70000,1.0e-3);
  plot(Ab2,Aa2,linestyle="-",color="black",linewidth=2)
end
legend()

ax1 = h1.axes[1]
ax1.set_xlabel(L"b", fontsize=lafs)
ax1.set_ylabel(L"a", fontsize=lafs)









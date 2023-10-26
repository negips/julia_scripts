#!/bin/julia

println("Testing Barkley's model")

using Roots
using PyPlot

include("Meinhardt_Nullclines.jl")        # For find_nullcline()
include("Barkley.jl")

close("all")

r = 0.8;

qnc_q,qnc_u,unc_q,unc_u = BarkleyNullClines(r);


n1 = size(qnc_q,2)
if n1>1
  lab1 = Vector{String}(undef,n1)
  lab1[1] = "∂q/∂t=0"
  for i in 2:n1
    lab1[i] = ""
  end
else
  lab1 = "∂q/∂t=0"
end  

n2 = size(unc_q,2)
if n2>1
  lab2 = Vector{String}(undef,n2)
  lab2[1] = "∂u/∂t=0"
  for i in 2:n2
    lab2[i] = ""
  end
else
  lab2 = "∂u/∂t=0"
end  

h1 = figure(num=1)
plot(qnc_u,qnc_q,linestyle="-",color="blue",linewidth=2, label=lab1)
plot(unc_u,unc_q,linestyle="-",color="red",linewidth=2,  label=lab2)
legend()


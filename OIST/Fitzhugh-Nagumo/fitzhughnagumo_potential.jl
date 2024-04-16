#!/bin/julia

println("Parametric Non-linearPotentials for the Fitzhugh-Nagumo Model")

using LinearAlgebra
using Random
using Printf
using PyPlot,PyCall
using Roots

include("FitzhughNagumo.jl")
include("../NullClines.jl")
include("$JULIACOMMON/MoveFigure.jl")


function cubic_potential_func(x,r)

  nr   = length(r)
  @assert nr == 3 "Only valid for cubic functions with 3 roots" 

  ϕ    = -( (1.0/4.0)*(x.^4) .- (1.0/3.0)*(x.^3)*(r[1] + r[2] + r[3]) .+ (1.0/2.0)*(x.^2)*(r[1]*r[2] + r[1]*r[3] + r[2]*r[3]) .- r[1]*r[2]*r[3]*x)

  return ϕ
end  
#---------------------------------------------------------------------- 

lafs  = 16
#a     = -0.50
#b     =  0.50

cm    = get_cmap("tab10");

close("all")

h1  = figure(num=1)
ax1 = h1.subplots()

cmapU = cm(0)
cmapV = cm(1)

Uu,Uv,Vu,Vv = FitzhughNagumoNullClines(a,b)

nullc_u = ax1.plot(Uu,Uv,linestyle="-",color=cmapU,linewidth=2)
nullc_v = ax1.plot(Vu,Vv,linestyle="-",color=cmapV,linewidth=2)

ax1.set_xlabel(L"u", fontsize=lafs)
ax1.set_ylabel(L"v", fontsize=lafs)

title = "a=$a; b=$b"
ax1.set_title(title, fontsize=lafs)
ax1.set_xlim(-3.0,3.0)
ax1.set_ylim(-3.0,4.0)
MoveFigure(h1,1250,500)


# vpar = LinRange(-1.99,1.99,9);
vpar = [-1.5; -1.0; 0.0; 1.0; 1.5];
nv   = length(vpar);
uroots = zeros(Float64,3,nv)

cm2    = get_cmap("tab20");

for i in 1:nv
  global uroots
  local v0  = vpar[i]
  ftilde(x) = FitzhughNagumo(x,v0,a,b)[1]
  rr        = find_zeros(ftilde,-3.0,3.0,verbose=false)
  r1        = minimum(rr)
  r2        = maximum(rr)
  # ax1.plot([r1; r2],[v0; v0],linestyle="--",color=cm2(i-1))
  if (length(rr)==3) 
    uroots[:,i] = rr
  end  

end  
  
include("plot_potential.jl")














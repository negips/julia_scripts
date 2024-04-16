#!/bin/julia

println("Calculate the phase flow of the Fitzhugh-Nagumo Model")

using LinearAlgebra
using Random
using Printf
using PyPlot,PyCall
using Roots

#include("Meinhardt.jl")
#include("Meinhardt_Nullclines.jl")
include("FitzhughNagumo.jl")
include("../NullClines.jl")
include("$JULIACOMMON/MoveFigure.jl")

lafs  = 16
a     = -0.50
b     =  0.50

cm    = get_cmap("tab10");

close("all")

h1 = figure(num=1)
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

#ax1.legend(fontsize=12)


#legend()


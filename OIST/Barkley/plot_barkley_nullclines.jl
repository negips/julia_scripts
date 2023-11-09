#!/bin/julia

println("Calculate the phase flow of the Barkley's Model for pipe flow")

#using LinearAlgebra
#using Random
#using Printf
#using PyPlot,PyCall
#using Roots

#include("Meinhardt.jl")
#include("Meinhardt_Nullclines.jl")
include("../NullClines.jl")

lafs = 16

cm    = get_cmap("tab10");

h1 = figure(num=1)
ax1 = h1.subplots()

cmapQ = cm(0)
cmapU = cm(1)

Qq,Qu,Uq,Uu = BarkleyNullClines(r)

nullc_q = ax1.plot(Qu,Qq,linestyle="-",color=cmapQ,linewidth=2)
nullc_u = ax1.plot(Uu,Uq,linestyle="-",color=cmapU,linewidth=2)

ax1.set_xlabel(L"u", fontsize=lafs)
ax1.set_ylabel(L"q", fontsize=lafs)

title = "R=$r"
ax1.set_title(title, fontsize=lafs)
ax1.set_ylim(-0.25,2.0)
ax1.set_xlim(0.5,4.0)
MoveFigure(h1,1250,500)

#ax1.legend(fontsize=12)


#legend()


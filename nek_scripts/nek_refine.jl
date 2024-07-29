#!/bin/julia

println("Refine 2 Elements in Nek (2D)")

using Statistics
using PyPlot

include("Nek_ElementRefine.jl")

ndim  = 2
nel   = 2
nc    = 2^ndim
nf    = 2*ndim

x0    = zeros(Float64,nc,nel)
y0    = zeros(Float64,nc,nel)

xmid  = zeros(Float64,nel)
ymid  = zeros(Float64,nel)

BC0   = fill("SYM",nf,nel)
Par0  = zeros(Int64,2,nf,nel) 

# 1st Element
e       = 1
x0[1,e] = -1.0
x0[2,e] =  0.0
x0[3,e] =  0.0
x0[4,e] = -1.0

y0[1,e] =  0.0
y0[2,e] =  0.0
y0[3,e] =  1.0
y0[4,e] =  1.0

f           = 1
BC0[f,e]    = "W  "
f           = 2
BC0[f,e]    = "E  "
Par0[1,f,e] = 2         # Connect to Element No.
Par0[2,f,e] = 4         # Connect on Face No.


# 2nd Element
e       = 2
x0[1,e] =  0.0
x0[2,e] =  1.0
x0[3,e] =  1.0
x0[4,e] =  0.0

y0[1,e] =  0.0
y0[2,e] =  0.0
y0[3,e] =  1.0
y0[4,e] =  1.0

BC0[1,e]  = "W  "
BC0[4,e]  = "E  "
f           = 1
BC0[f,e]    = "W  "
f           = 4
BC0[f,e]    = "E  "
Par0[1,f,e] = 1         # Connect to Element No.
Par0[2,f,e] = 2         # Connect on Face No.


xmid = mean(x0,dims=1)
ymid = mean(y0,dims=1)

close("all")

plot(x0,y0,linestyle="-",linewidth=2,marker="o")

# Refine Bottom Right
x     = zeros(Float64,4,6)
y     = zeros(Float64,4,6)
BC    = fill("E  ",nf,6)
Par   = zeros(Int64,2,nf,6)

x[:,1:3],y[:,1:3],BC[:,1:3],Par[:,:,1:3] = Nek_ElementRefine(x0[:,1],y0[:,1],BC0[:,1],Par0[:,:,1],2)
x[:,4:6],y[:,4:6],BC[:,4:6],Par[:,:,4:6] = Nek_ElementRefine(x0[:,2],y0[:,2],BC0[:,2],Par0[:,:,2],1)

plot(x,y,linestyle="-.")

println("Done.")









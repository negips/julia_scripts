#!/bin/julia

println("Testing Writing to file")

include("../../SEMla/Module_SEMla/Module_Nek/Nek.jl")

ver   = 19.0
ndim  = 2
nelg  = 1
nelv  = 1
ifre2 = false
npar  = 5

# Open File
io    = open("test.rea", "w")

Nek.write_nekhdr(io,ver,ndim)

pars = zeros(Float64,200)
P,C  = Nek.param_descriptions()
Nek.write_params(io,pars[1:npar],P[1:npar],C[1:npar])
Nek.write_ps_hdr(io,0)
Nek.write_logicalswitches_hdr(io,0)
Nek.write_prenek_hdr(io)
Nek.write_mesh_hdr(io,ndim,nelg,nelv,ifre2)

XC    = rand(Float64,4,2)
YC    = rand(Float64,4,2)
Nek.write_mesh_coords(io,ndim,XC,YC,XC)
ncurve = 0
Nek.write_curve_hdr(io,ncurve)
Nek.write_fluidBC_hdr(io)

close(io)

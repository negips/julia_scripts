#!/bin/julia

println("Testing Writing to file")

include("../../SEMla/Module_SEMla/Module_Nek/Nek.jl")

ver   = 19.0
ndim  = 2
nelg  = 2
nelv  = 2
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
ZC    = zeros(Float64,4,2)

BC0   = fill("SYM",4,2)
Par0  = zeros(Int64,5,4,2) 

Nek.write_mesh_coords(io,ndim,XC,YC,XC)
ncurve = 0
Nek.write_curve_hdr(io,ncurve)
Nek.write_fluidBC_hdr(io)
Nek.write_fluidBC(io,ndim,nelg,BC0,Par0)

Nek.write_zerothermalBC_hdr(io)

nrst     = 0
rstfiles = Vector{String}(undef,nrst)
rstopts  = Vector{String}(undef,nrst)
Nek.write_restart(io,nrst,rstfiles,rstopts)

nic     = 0
icfiles = Vector{String}(undef,nic)
Nek.write_ic(io,nic,icfiles)

ndf     = 0
driveforce = Vector{String}(undef,ndf)
Nek.write_driveforce(io,ndf,driveforce)

nvp         = 1
npackets    = 0
datapacket  = Vector{String}(undef,0)
Nek.write_varprop(io,nvp,npackets,datapacket)

nhist = 0
history = Vector{String}(undef,nhist)
Nek.write_hist(io,nhist,history)

niospec     = 6
OutFlds     = Vector{String}(undef,niospec)
OutFlds[1]  = "COORDINATES"
OutFlds[2]  = "VELOCITY"
OutFlds[3]  = "PRESSURE"
OutFlds[4]  = "TEMPERATURE"
OutFlds[5]  = "TEMPERATURE GRADIENT"
OutFlds[6]  = "PASSIVE SCALARS"
Specs       = Vector{String}(undef,niospec)
Specs[1]    = "T"
Specs[2]    = "T"
Specs[3]    = "T"
Specs[4]    = "F"
Specs[5]    = "F"
Specs[6]    = "0"

Nek.write_outputspec(io,niospec,OutFlds,Specs)

nobj        = 4
Objects     = Vector{String}(undef,nobj)
Objects[1]  = "Surface  Objects"
Objects[2]  = "Volume   Objects"
Objects[3]  = "Edge     Objects"
Objects[4]  = "Point    Objects"
ObjSpecs    = zeros(Int64,nobj)

Nek.write_objectspec(io,nobj,Objects,ObjSpecs)

close(io)












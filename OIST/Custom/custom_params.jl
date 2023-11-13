#!/bin/julia

# Setting the parameters here

U                 = prec(0)         # Convection
γ                 = prec(1)         # Diffusion (Generic) 
γa                = prec(0.2)       # Diffusion (activator)
γb                = prec(1.0)       # Diffusion (inhibitor)
γall              = [γb γa]
σa                =  0.0             # Activator Noise Strength
σb                =  0.0             # Inhibitor Noise Strength
σall              = [σb σa]

ifsparse          = true
ifperiodic        = true
ifglobal          = true

# Initialization
ngauss            = 1
x0gauss           = [20.0 80.0]           #xe*rand(ngauss)
ampgauss          = ones(Float64,ngauss)  #rand(ngauss)
x0                = x0gauss[1]      # Gaussian Center
σg                = prec(2)         # Gaussian Std. Deviation
ampA0             = prec(0.0)
ampB0             = prec(0)

A0Off             = prec(1)         # Homogeneous state value
B0Off             = prec(0)         # Homogeneous state value

σai               = prec(3.3)       # Initial Activator Noise
σbi               = prec(0.0)         # Initial Inhibitor Noise Strength


# Simulation
dt                = prec(0.01)
nsteps            = 12000
nstep_switch1     = 3000 
nstep_switch2     = 12100 

verbosestep       = 100
plotupd           = 50
ifphplot          = true
surf_save         = 20
nsurf_save        = floor(Int,nsteps/surf_save)+1

ifpldyn           = false      # Plot dynamic phase


# Saving Params
# Standard diffusion params used:
# γa = 0.2
# γb = 0.01

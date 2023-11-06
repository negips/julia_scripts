#!/bin/julia

# Setting the parameters here


U                 = prec(0)         # Convection
γ                 = prec(1)         # Diffusion (Generic) 
γa                = prec(0.2)       # Diffusion (activator)
γb                = prec(0.01)      # Diffusion (inhibitor)
γall              = [γb γa]
σa                = 10.0            # Activator Noise Strength
σb                = 5.0             # Inhibitor Noise Strength
σall              = [σb σa]

ifsparse          = true
ifperiodic        = true
ifglobal          = true

# Initialization
x0                = prec(100)       # Gaussian Center
σg                = prec(2)         # Gaussian Std. Deviation
ampA0             = prec(4.0)
ampB0             = prec(0)

A0Off             = prec(0)         # Homogeneous state value
B0Off             = prec(0)         # Homogeneous state value

σai               = prec(0.5)       # Initial Activator Noise
σbi               = prec(0)         # Initial Inhibitor Noise Strength


# Simulation
dt                = prec(0.01)
nsteps            = 20000
nstep_switch1     = 3000 
nstep_switch2     = 3100 

verbosestep       = 100
plotupd           = -1
ifphplot          = true
surf_save         = 100
nsurf_save        = floor(Int,nsteps/surf_save)+1

ifpldyn           = true      # Plot dynamic phase


# Saving Params
# Standard diffusion params used:
# γa = 0.2
# γb = 0.01

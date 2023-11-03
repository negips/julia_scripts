#!/bin/julia

# Setting the parameters here


U                 = prec(0)         # Convection
γ                 = prec(1)         # Diffusion (Generic) 
γa                = prec(0.2)       # Diffusion (activator)
γb                = prec(0.2)       # Diffusion (inhibitor)
γall              = [γb γa]
σa                = 0.0             # Activator Noise Strength
σb                = 0.0             # Inhibitor Noise Strength
σall              = [σb σa]

ifsparse          = true
ifperiodic        = true
ifglobal          = true

# Initialization
x0                = prec(25)       # Gaussian Center
σg                = prec(2)         # Gaussian Std. Deviation
ampA0             = prec(2)
ampB0             = prec(0)

A0Off             = prec(0)         # Homogeneous state value
B0Off             = prec(0)         # Homogeneous state value

σai               = 10.0            # Initial Activator Noise
σbi               = 10.0            # Initial Inhibitor Noise Strength


# Simulation
dt                = prec(0.01)
nsteps            = 15000
nstep_switch      = 7000 

verbosestep       = 1000
plotupd           = 50
ifphplot          = true
surf_save         = 100
nsurf_save        = floor(Int,nsteps/surf_save)+1




# Saving Params
# Standard diffusion params used:
# γa = 0.2
# γb = 0.01

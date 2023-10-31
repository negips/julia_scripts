#!/bin/julia

# Setting the parameters here


U           = prec(0)         # Convection
γ           = prec(1)         # Diffusion (Generic) 
γa          = prec(0.2)       # Diffusion (activator)
γb          = prec(0.01)       # Diffusion (inhibitor)
γall        = [γb γa]

ifsparse    = true
ifperiodic  = true
ifglobal    = true

# Initialization
x0          = prec(100)       # Gaussian Center
σ           = prec(2)         # Gaussian Std. Deviation
ampA0       = prec(1)
ampB0       = prec(0.0)

A0Off       = prec(0)         # Homogeneous state value
B0Off       = prec(0)         # Homogeneous state value


# Simulation
dt          = prec(0.01);
nsteps      = 100000;

plotupd     = 100;
ifphplot    = true

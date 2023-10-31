#!/bin/julia

# Setting the parameters here


U           = prec(0)         # Convection
γ           = prec(1)         # Diffusion (Generic) 
γa          = prec(0.015)     # Diffusion (activator)
γb          = prec(0.00)      # Diffusion (inhibitor)
γall        = [γa γb]

ifsparse    = true
ifperiodic  = true
ifglobal    = true

R           = prec(0.01)       # Meinhardt model parameter

# Initialization
x0          = prec(100)       # Gaussian Center
σ           = prec(1)         # Gaussian Std. Deviation
ampA0       = prec(3)
ampB0       = prec(1)

A0Off       = prec(0)         # Homogeneous state value
B0Off       = prec(1)         # Homogeneous state value


# Simulation
dt          = prec(0.01);
nsteps      = 100000;

plotupd     = 50;

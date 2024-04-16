#!/bin/julia

# Setting the parameters here

U           = prec(1)         # Convection Coeff.
γ           = prec(1)         # Diffusion Coeff. (Generic) 
μ           = prec(1)         # Source Coeff.

# Diffusion parameters
γu          = prec(0.01)      # Diffusion (u)
γq          = prec(0.50)      # Diffusion (q)
γall        = [γu γq]

# Convection parameters
ζu          = prec(0)
ζq          = prec(0.8)
ζall        = [ζu ζq]

# Artificially added. Speed up factor for convting field (u).
# 1.0 => Normal speed
cfu         = prec(1.0)
cfq         = prec(1.0)
cfall       = [cfu cfq]


ifsparse    = true
ifperiodic  = true
ifglobal    = true

r           = prec(0.6)       # Barkley Reynolds number

# Initialization
x0          = prec(20)       # Gaussian Center
σ           = prec(2)         # Gaussian Std. Deviation
ampU0       = prec(0)
ampQ0       = prec(1)

U0Off       = prec(2)         # Homogeneous state value
Q0Off       = prec(0)         # Homogeneous state value

# Simulation
dt          = prec(0.003);
nsteps      = 50000;

plotupd     = 100
ifphplot    = true

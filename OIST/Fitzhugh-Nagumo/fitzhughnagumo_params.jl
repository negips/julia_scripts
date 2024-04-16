#!/bin/julia

# Setting the parameters here

U           = prec(1)         # Convection Coeff.
γ           = prec(1)         # Diffusion Coeff. (Generic) 
μ           = prec(1)         # Source Coeff.

# Diffusion parameters
γu          = prec(0.10)       # Diffusion (u)
γv          = prec(0.00)      # Diffusion (v)
γall        = [γu γv]

# Convection parameters
ζu          = prec(0.0)
ζv          = prec(0.0)
ζall        = [ζu ζv]

# Artificially added. Speed up factor for convting field (u).
# 1.0 => Normal speed
cfu         = prec(1.0)
cfv         = prec(1.0)
cfall       = [cfu cfv]


ifsparse    = true
ifperiodic  = true
ifglobal    = true

a           = prec(-0.5)       #
b           = prec(0.5)        #
ϵ           = 0.1


# Initialization
x0          = prec(25)       # Gaussian Center
σ           = prec(1)         # Gaussian Std. Deviation
ampU0       = prec(2)
ampV0       = prec(0)

U0Off       = prec(-1.2)         # Homogeneous state value
V0Off       = prec(-2)         # Homogeneous state value

# Simulation
dt          = prec(0.003);
nsteps      = 20000;

plotupd     = 100
ifphplot    = true

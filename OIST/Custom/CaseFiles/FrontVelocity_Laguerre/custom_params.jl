#!/bin/julia

# Setting the parameters here

U                 = prec(1)         # Convection
γ                 = prec(1)         # Diffusion (Generic) 
γa                = prec(0.2)       # Diffusion (activator)
γb                = prec(0.01)      # Diffusion (inhibitor)
γζ                = prec(0.001)        # Diffusion for auxillary variable
γall              = [γb γa γζ]
σa                =  0.0e-2          # Activator Noise Strength
σb                =  0.0             # Inhibitor Noise Strength
σζ                =  0.0             # Aux. Noise Strength
σall              = [σb σa σζ]

ifsparse          = false
ifperiodic        = false
ifglobal          = true

nflds             = 2               # No of fields

# Initialization
ngauss            = 1
x0gauss           = [80.0]           #xe*rand(ngauss)
ngauss            = length(x0gauss)
ampgauss          = ones(Float64,ngauss)  #rand(ngauss)
x0                = x0gauss[1]      # Gaussian Center
σg                = prec(2)         # Gaussian Std. Deviation
ampA0             = prec(1.0)
ampB0             = prec(0.0)
ampζ0             = prec(0.0)

Amp0              = zeros(prec,nflds) 
Amp0              = [ampB0; ampA0; ampζ0]

A0Off             = prec(0.0)         # Homogeneous state value
B0Off             = prec(0.0)         # Homogeneous state value
ζ0Off             = prec(0.0)         # Homogeneous state value

Off0              = zeros(prec,nflds)
Off0              = [B0Off; A0Off; ζ0Off]

σai               = prec(0.0)         # Initial Activator Noise
σbi               = prec(0.0)         # Initial Inhibitor Noise Strength
σζi               = prec(0.0)         # Initial Aux. Strength
σ0i               = [σbi; σai; σζi]


# Simulation
dt                = prec(0.01)       # Time step
nsteps            = 1000             # No of steps
nstep_switch1     = 3000             # Switch functions 1
nstep_switch2     = 3400             # Switch again

Ω                 = 0.02             # Time scale of parameter oscillation

verbosestep       = 100
plotupd           = 20
surf_save         = 20
nsurf_save        = floor(Int,nsteps/surf_save)+1


iffldplot         = true      # Plot fields
ifphplot          = true      # Plot Phase A-B
ifdynplot         = false     # Plot dynamic phase (ζ)
initplot          = true      # Plot initial conditions
ifdynnull         = true      # Plot null-clines dynamically
ifsaveframe       = false


ifplot            = iffldplot || ifphplot || ifdynplot 
plotfldi          = fill(true,nflds)
plotfldi[1]       = true


# Saving Params
# Standard diffusion params used:
# γa = 0.2
# γb = 0.01




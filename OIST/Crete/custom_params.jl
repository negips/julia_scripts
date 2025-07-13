#!/bin/julia

# Setting the parameters here
# FitzHugh-Nagumo Params
a                 =  5.0            # Third root
b                 = -4.0            #
R                 =  1.0            # Strength of g
m                 =  1.50           # Slope of null-cline
ϵ                 =  1.0E-1
η                 =  1.0E+0

U                 = 0               # Convection
γ                 = 1               # Diffusion (Generic) 
γ1                = 1.0*ϵ^2         # Diffusion (v)
γ2                = 1.0*ϵ           # Diffusion (u)
γ3                = 0.00            # Diffusion auxillary variable
γall              = [γ1 γ2 γ3]
σ1                =  0.0e-2         # Noise Strength v
σ2                =  0.0            # Noise Strength u
σ3                =  0.0            # Aux. Noise Strength
σall              = [σ1 σ2 σ3]

ifsparse          = true
ifperiodic        = false
ifglobal          = true

nflds             = 2               # No of fields

# Initialization
ngauss            = 0
x0gauss           = [25.0]           # xe*rand(ngauss)
ngauss            = length(x0gauss)
ampgauss          = ones(Float64,ngauss)  #rand(ngauss)
x0                = x0gauss[1]      # Gaussian Center
σg                = 2               # Gaussian Std. Deviation
ampA0             = 2.0
ampB0             = 0.0
ampζ0             = 0.0

Amp0              = zeros(Float64,nflds) 
Amp0              = [ampB0; ampA0; ampζ0]

Normk0            = [8.0]           # Normalized Wavenumbers
nk0               = length(Normk0)
Ak0               = [0.0; 0.25]
Bk0               = [0.0; 0.25]
ampk0             = [Ak0 Bk0]


A0Off             = 0.0             # Homogeneous state value
B0Off             = 0.0             # Homogeneous state value
ζ0Off             = 0.0             # Homogeneous state value

Off0              = zeros(Float64,nflds)
Off0              = [B0Off; A0Off; ζ0Off]

σai               = 0.0             # Initial Activator Noise
σbi               = 0.0             # Initial Inhibitor Noise Strength
σζi               = 0.0             # Initial Aux. Strength
σ0i               = [σbi; σai; σζi]


# Simulation
dt                = 2.0e-3           # Time step
nsteps            = 5000            # No of steps

verbosestep       = 100
plotupd           = 10
surf_save         = 10
nsurf_save        = floor(Int,nsteps/surf_save)+1


iffldplot         = true            # Plot fields
ifphplot          = true            # Plot Phase A-B
ifdynplot         = false           # Plot dynamic phase (λ)
initplot          = true            # Plot initial conditions
ifdynnull         = false           # Plot null-clines dynamically
ifsaveframe       = true            # save individual frames
ifsavext          = true            # Save X-T plot


ifplot            = iffldplot || ifphplot || ifdynplot 
plotfldi          = fill(true,nflds)
plotfldi[1]       = true

Aeq               = 0.50
Asen              = 2.50


# Saving Params
# Standard diffusion params used:
# γa = 0.2
# γb = 0.01




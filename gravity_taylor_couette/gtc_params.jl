#!/usr/bin/julia

println("Parameters for the Gravity Taylor Couette Problem")
println("Marcado et al. (2015) Gravity Wave instability in aturbulent free-surface Taylor-Couette flow: experiments and comparison with an amplitude equation with additive noise, New J. Phys. 17, 013039 ")

using LinearAlgebra

# Dimensional values follow notation of paper
# Units and descriptions in comments.

println(" ")

ri = 0.020                    # (m) Inner Radius
ro = 0.166                    # (m) Outter Radius

Hc = 0.600                    # (m) Total Cylinder Height

g  = 9.81                     # (m/s2) Gravitational Acceleration
L  = ro-ri                    # (m) Cylinder Gap

h0 = 1.23*L                   # (m) Height of the water column

k  = 1.0                      # (1/m) Azimuthal Wavenumber
#ω  = sqrt(g*k*tanh(g*h0))    # (1/s) Fundamental Mode Angular Frequency
f0 = 1.0/0.625                # (1/s) Fundamental Mode Frequency
τ0 = 1.0/f0                   # (s) Time Period of the Fundamental Mode
ω  = 2.0*π*f0                 # (1/s) Fundamental Mode Angular Frequency

Ωm = 10.4*ω                   # (1/s) Angular Rotation Rate of Inner Cylinder
fm = Ωm/(2.0*π)               # (1/s) Rotation Rate of Inner Cylinder
U  = Ωm*ri                    # (m/s) Azimuthal Velocity of Inner Cylinder

# # 4° C Values
# ρ  = 1000.0                   # (kg/m3) Density of Water at 4° C
# ν  = 1.5707e-6                # (m2/s) Kinematic Viscosity of Water at 4° C
# μ  = ρ*ν                      # (kg/m-s) Dynamic Viscosity of Water at 4° C
# σ  = 0.0751                   # (kg/s2) Surface Tension of Water at 4° C

# 20° C Values for Water
ρ  = 998.19                   # (kg/m3) Density of Water at 20° C
ν  = 1.0023e-6                # (m2/s) Kinematic Viscosity of Water at 20° C
μ  = ρ*ν                      # (kg/m-s) Dynamic Viscosity of Water at 4° C
σ  = 0.0728                   # (kg/s2) Surface Tension of Water at 4° C

# 20° C Values for Air
ρa  = 1.204                    # (kg/m3) Density of Air at 20° C
νa  = 15.06e-6                 # (m2/s) Kinematic Viscosity of Air at 20° C
μa  = ρ*ν                      # (kg/m-s) Dynamic Viscosity of Air at 4° C
σa  = 0.0                      # (kg/s2) Surface Tension of Air at 4° C

# Non-Dimensional Parameters
Re_w    = U*L/ν                           # Reynolds Number
Fr_w    = (U*U)/(g*h0)                    # Froude Number
Ca_w    = ρ*ν*U/σ                         # Capilary Number
h0byL   = h0/L                            # Water Height
η0      = ri/ro                           # Radii ratio
Ta_w    = ((Ωm^2)*ri*(ro-ri)^3)/ν^2       # Taylor number

HcbyL   = Hc/L                            # Total Column Height

# Non-Dimensional values for Air
Re_a    = U*L/νa                 # Reynolds Number
Ta_a    = ((Ωm^2)*ri*(ro-ri)^3)/νa^2       # Taylor number


# Non-dimensionalization

Re    = 3.0e5                 # Reynolds Number
Fr    = 1.0                   # Froude Number
Ca    = Ca_w                  # Capilary Number
h0byL = 1.23                  # Water Height
η0    = ri/ro                 # Radii ratio

U_nek  = 1.0
L_nek  = 1.0
ρ_nek  = 1.0
ri_nek = 1.0

h0_nek = h0byL*L_nek
ν_nek  = U_nek*L_nek/Re
g_nek  = U_nek*U_nek/h0_nek
σ_nek  = ρ_nek*ν_nek*U_nek/Ca
ro_nek = ri_nek/η0
hc_nek = HcbyL*L_nek







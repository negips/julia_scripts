#!/bin/julia

println("Building a null cline based on given Points/gradients by Minimizing the Lagrangian (Error)")


using LinearAlgebra
using Roots
using PyPlot
using PolynomialBases

const SRC = "/home/prabal/workstation/git/julia/OIST/Custom"

include("$SRC/BuildTimeDerivatives.jl")
include("$SRC/NullClines.jl")
include("$SRC/NullClineFcn.jl")
include("$SRC/ElementOf.jl")
include("$SRC/NullClineParams.jl")
include("$JULIACOMMON/MoveFigure.jl")

close("all")

lafs = 16

#include("select_nullclines.jl")
sets              = [2]

cm                = get_cmap("tab10")

set               = sets[1]
pars              = GetNullClineParams(set) 

h1                = figure(num=1)
MoveFigure(h1,1250,830)

# Container to hold plot handles
PlotContainers    = Array{Any}(undef,20)

Basis             = LobattoLegendre(1)


# Activator
#---------------------------------------------------------------------- 

f(x,y)            = FXY(x,y,pars.fc0,pars.fcx,pars.fcy)
if f(0.0,100.0)>0
  pars.fc0        = -pars.fc0
  pars.fcx        = -pars.fcx
  pars.fcy        = -pars.fcy
end  

# If we want +λ to be stabilizing or destabilizing
stabilizing = true
# Translated
if stabilizing 
  ϕfd   = 120 #150.0
else
  ϕfd   = -80.0
end  
println("F(x,y) Translated with Slope: $ϕfd Degrees")
ϕf      = ϕfd*π/180.0

# Plot the null-cline
λ0          = 0.0
dλ          = 0.35
#λvalues     = [λ0-dλ; λ0; λ0+dλ]
nλ          = 10
r_values    = LinRange(-1.0,1.0,nλ)
λmin        = λ0 - dλ
λmax        = λ0 + dλ
#λvalues     = LinRange(λ0-dλ,λ0+dλ,10)
Nullfx      = Vector{Vector{Float64}}(undef,nλ)
Nullfy      = Vector{Vector{Float64}}(undef,nλ)

e = 0
for r in r_values
  λ                   = map_from_canonical(r,λmin,λmax,Basis)
  # local f2(x,y)       = RotFXY(x,y,θ0,pars.fc0,pars.fcx,pars.fcy)
  local f2(x,y)       = TransFXY(x,y,λ,ϕf,pars.fc0,pars.fcx,pars.fcy)
  # Plot the null-cline
  local xi            = -20.0
  local yr0           = -10.0
  local yr1           =  50.0
  local dτ            = 1.0e-3
  local nsteps        = 200000
  local f0x,f0y       = NullClines(f2,xi,yr0,yr1,nsteps,dτ)
  global e           += 1
  Nullfx[e]           = f0x[:]
  Nullfy[e]           = f0y[:]
#  PlotContainers[plc] = ax1.plot(f0x,f0y,linestyle="-",label="λ=$λ; ϕ=$ϕfd")
end  


# Inhibitor
#---------------------------------------------------------------------- 
g(x,y)            = FXY(x,y,pars.gc0,pars.gcx,pars.gcy)
if g(100.0,0.0)>0
  pars.gc0        = -pars.gc0
  pars.gcx        = -pars.gcx
  pars.gcy        = -pars.gcy
  g(x,y)          = FXY(x,y,pars.gc0,pars.gcx,pars.gcy)
end  
stabilizing = true
# Translated
if stabilizing
  ϕgd             = -30.0
else
  ϕgd             = 150.0
end  
#println("G(x,y) Translated with Slope: $ϕgd Degrees")
ϕg                = ϕgd*π/180.0

θ0          =  0.0
dθ          =  5.0
nθ          = 10
s_values    = LinRange(-1.0,1.0,nθ)
θmin        = θ0 - dθ
θmax        = θ0 + dθ

Nullgx  = Vector{Vector{Float64}}(undef,nθ)
Nullgy  = Vector{Vector{Float64}}(undef,nθ)

e = 0
for s in s_values
  local θ             = (pi/180.0)*map_from_canonical(s,θmin,θmax,Basis)
  local g2(x,y)       = RotFXY(x,y,θ,pars.gc0,pars.gcx,pars.gcy)
  local xi            = -10.0
  local yr0           = -30.0
  local yr1           =  10.0
  local dτ            = 1.0e-3
  local nsteps        = 200000
  local g20x,g20y     = NullClines(g2,xi,yr0,yr1,nsteps,dτ)
  global e           += 1
  Nullgx[e]           = g20x[:]
  Nullgy[e]           = g20y[:]
  # PlotContainers[plc] = ax1.plot(g20x,g20y,linestyle="--",label="θ=$λ; ϕ=$ϕgd")
  # legend()
end  


l1 = minimum(length.(Nullfx))

Xf  = zeros(Float64,l1,nλ)
Yf  = zeros(Float64,l1,nλ)
Zf  = zeros(Float64,l1,nλ)

for i in 1:nλ
  Xf[:,i] = Nullfx[i][1:l1]
  Yf[:,i] = Nullfy[i][1:l1]
  z       = r_values[i]
  local λ = map_from_canonical(z,λmin,λmax,Basis)
  @views fill!(Zf[:,i],z)
end

cm2 = get_cmap("Blues")
surf(Xf,Yf,Zf,cmap=cm2,edgecolor="none")

l2 = minimum(length.(Nullgx))

Xg  = zeros(Float64,l2,nθ)
Yg  = zeros(Float64,l2,nθ)
Zg  = zeros(Float64,l2,nθ)

for i in 1:nθ
  Xg[:,i] = Nullgx[i][1:l2]
  Yg[:,i] = Nullgy[i][1:l2]
  z       = s_values[i]
  local θ = map_from_canonical(z,θmin,θmax,Basis)
  @views fill!(Zg[:,i],z)
end

cm3 = get_cmap("OrRd") 
surf(Xg,Yg,Zg,cmap=cm3,edgecolor="none")
ax = h1.gca()
ax.roll = 90;
ax.elev = 90;
ax.azim = 0; 
draw()
#ax.set_xlim(-2.0,7.0)
#ax.set_ylim(-5.0,10.0)

#for e in 1:nλ
#  plot3D(Nullfx[e],Nullfy[e],r_values[e])
#end
#
#for e in 1:nθ
#  plot3D(Nullgx[e],Nullgy[e],s_values[e])
#end
#
#MoveFigure(h1,1250,830)
#ax = h1.gca()









println("Main interface for Time Stepper")
println("Using Implicitly Restarted Arnoldi")

using PyPlot,PyCall
using LinearAlgebra
using IterativeSolvers
using SpecialFunctions
using Roots
using Random
# using GenericLinearAlgebra          # For eigvals for BigFloat
# using DoubleFloats
using Printf
using JLD2


include("ArnUpd.jl")
include("ArnIRst.jl")
include("ExplicitShiftedQR.jl")
include("BulgeChase.jl")
include("IRAM.jl")
#include("RK4.jl")
include("$JULIACOMMON/RK4.jl")

close("all")

# Ifglobal
ifglobal = true

rcParams = PyPlot.PyDict(PyPlot.matplotlib."rcParams")

# Analytical Eigenvalues
ω1    = find_zero(airyai,(-3.0,-0.0))
ω2    = find_zero(airyai,(-5.0,-3.0))
ω3    = find_zero(airyai,(-6.0,-5.0))
ω4    = find_zero(airyai,(-7.0,-6.0))
ω5    = find_zero(airyai,(-8.0,-7.0))
ω6    = find_zero(airyai,(-9.5,-8.0))
ω7    = find_zero(airyai,(-10.5,-9.5))
ω8    = find_zero(airyai,(-11.8,-10.5))
ω9    = find_zero(airyai,(-12.0,-11.8))
ω10   = find_zero(airyai,(-12.9,-12.0))
ω11   = find_zero(airyai,(-13.8,-12.9))
ω12   = find_zero(airyai,(-14.8,-13.8))
ω13   = find_zero(airyai,(-15.8,-14.8))
ω14   = find_zero(airyai,(-16.8,-15.8))
ω15   = find_zero(airyai,(-17.5,-16.8))

ω     = [ω1, ω2, ω3, ω4, ω5, ω6, ω7, ω8, ω9, ω10, ω11, ω12, ω13, ω14, ω15]

Ω0    = im*1.0
U     = 1.0                                           # -δ_1
ϕ     = -1.0*π/4.0
R     = 1.0
γ     = R*exp(im*ϕ)                                   # δ_4
μx    = -U/8.0                                        # δ_3
μ0    = Ω0 + (U^2)/(4.0*γ) - ((γ*μx*μx)^(1.0/3.0))*ω1 # δ_2


δγ    = 0.0*(exp(im*0))
γ     = γ + δγ
δμx   = -0.00*(U/8)
μx    = μx + δμx

ifconj = false
if (ifconj)
  U    = conj(U)        # -δ_1
  γ    = conj(γ)        # δ_4
  μ0   = conj(μ0)       # δ_2
  μx   = conj(μx)       # δ_3
end  

Ω  = (μ0 .- U*U/(4.0*γ) .+ (γ*μx*μx)^(1.0/3.0)*ω)

# Plots
#-------------------------------------------------- 
rcParams["markers.fillstyle"] = "none"
hλ = figure(num=1,figsize=[8.,6.]);
ax1 = gca()
pΛ = ax1.plot(imag.(Ω),real.(Ω),linestyle="none",marker="o",markersize=8)
Ωi_max = maximum(abs.(imag.(Ω)))
Ωr_min = minimum(real.(Ω))
Ωr_max = maximum(real.(Ω))

ax1.set_xlim((-2*Ωi_max,2*Ωi_max))
ax1.set_ylim((1.1*Ωr_min,1.0))
rcParams["markers.fillstyle"] = "full"
# Include the function files


include("sem_main.jl")
xg          = QT*(vimult.*Geom.xm1[:])
OPg[1,:]   .= 0.

Bgi   = 1.0./Bg
ev2   = eigvals(Matrix(OPg),diagm(Bg))
pl2   = ax1.plot(imag.(ev2),real.(ev2),linestyle="none",marker="s",markersize=8)







println("Done.")








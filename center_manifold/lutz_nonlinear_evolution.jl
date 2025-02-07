println("Non-Linear evolution for Ginzburg Landau equations")

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
include("NLGinzburgLandau.jl")
include("OP_RK4.jl")

close("all")

# Ifglobal
ifglobal = true

rcParams = PyPlot.PyDict(PyPlot.matplotlib."rcParams")

# Analytical Eigenvalues
ω1 = find_zero(airyai,(-3.0,-0.0))
ω2 = find_zero(airyai,(-5.0,-3.0))
ω3 = find_zero(airyai,(-6.0,-5.0))
ω4 = find_zero(airyai,(-7.0,-6.0))
ω5 = find_zero(airyai,(-8.0,-7.0))
ω6 = find_zero(airyai,(-9.5,-8.0))
ω7 = find_zero(airyai,(-10.5,-9.5))
ω8 = find_zero(airyai,(-11.8,-10.5))
ω9 = find_zero(airyai,(-12.0,-11.8))
ω10 = find_zero(airyai,(-12.9,-12.0))
ω11 = find_zero(airyai,(-13.8,-12.9))
ω12 = find_zero(airyai,(-14.8,-13.8))
ω13 = find_zero(airyai,(-15.8,-14.8))
ω14 = find_zero(airyai,(-16.8,-15.8))
ω15 = find_zero(airyai,(-17.5,-16.8))

ω  = [ω1, ω2, ω3, ω4, ω5, ω6, ω7, ω8, ω9, ω10, ω11, ω12, ω13, ω14, ω15]
#Ω  = im*(U*U/8.0 .- U*U/(4.0*γ) .+ γ^(1.0/3.0)*(U^(4.0/3.0))/(160.0^(2.0/3.0))*ω)

# Ginzburg Landau Parameters
Ω0    = im*1.0
U     = 1.0
ϕ     = -π/4.0
R     = 1.0
γ     = R*exp(im*ϕ)
μx    = U/8.0 
μ0    = Ω0 + (U^2)/(4.0*γ) - ((γ*μx*μx)^(1.0/3.0))*ω1
δ5    = -0.1 + 0.1im

Ω  = (μ0 .- U*U/(4.0*γ) .+ (γ*μx*μx)^(1.0/3.0)*ω)

# Include the function files
include("sem_main.jl")

rng   = MersenneTwister(1235)

xg    = QT*(vimult.*Geom.xm1[:])
vt    = VT

v     = randn(vt,ndof)*0.001;

ifplot      = true 
verbose     = true
nsteps      = 500000
ifsave      = false
plotstep    = 1000
verbosestep = 1000
#if (ifadjoint)
#  Ω = conj.(Ω)
#end  
if (ifplot)
  hv  = figure(num=2,figsize=[8.,6.]);
  ax2 = gca()
end

cm    = get_cmap("tab10");
rgba0 = cm(0) 
rgba1 = cm(1) 
rgba2 = cm(2) 

dt    = prec(0.0001)

bc    = zeros(vt,1,ndof);
NGL(x)= NLGinzburgLandau(OPg,ones(vt,ndof),x,δ5,vt(0),vt(0),true,false)  

t = zro*dt        # Time
i = 0             # Istep

maxouter_it = 150
major_it    = 1

# Start iterations
ifdirect = true
println("Starting Iterations")

ifconv = false
while (~ifconv)
  global v
  global t, i
  global plr,pli
  global OPg
  global hλ, ax1
  global ax2

  local β


  β = one
  i = i + 1

  t = t + dt;

## Build the operator only the first time
#  if (i==1)
##   Direct Operator BCs 
#    OPg[1,:] = bc
#    OPg[1,1] = one + im*zro        # Change operator for BC
#  end  

# Apply BC
  v[1]    = zro + im*zro
  v       = OP_RK4!(NGL,v,dt)

# No Arnoldi iteration      
  if verbose && mod(i,verbosestep)==0
    println("Istep=$i, Time=$t")
  end

# Plot the field  
  if (ifplot && mod(i,plotstep)==0)
    if (i>plotstep) 
      for lo in ax2.get_lines()
        lo.remove()
      end  
    end  
   
    pv1 = ax2.plot(xg,real.(v),linestyle="-",color=rgba0)

    vmin = 2.0*minimum(real.(v))
    vmax = 2.0*maximum(real.(v))
    dv   = abs(vmax-vmin)
#    ax2.set_ylim((vmin,vmax))
    ax2.set_ylim((-dv,dv))
   
  end 
  
  if i==nsteps
    break
  end  

end       # while ... 


println("Done.")








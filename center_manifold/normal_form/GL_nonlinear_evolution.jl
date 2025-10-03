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
using Peaks
using Statistics

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

lafs        = 16
lgfs        = 12

# Ifglobal
ifglobal    = true

rcParams    = PyPlot.PyDict(PyPlot.matplotlib."rcParams")

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
#Ω  = im*(U*U/8.0 .- U*U/(4.0*γ) .+ γ^(1.0/3.0)*(U^(4.0/3.0))/(160.0^(2.0/3.0))*ω)

# Ginzburg Landau Parameters
Ω0    = 0.0
U     = 1.0
ϕ     = -0.0*π/4.0
R     = 1.0

γ     = R*exp(im*ϕ)

μx    = U/8.0 
μ0    = Ω0 + (U^2)/(4.0*γ) - ((γ*μx*μx)^(1.0/3.0))*ω1
#δ5    = (-1.0 + 1.0im)*0.1
δ5    = (-1.0 + 0.0im)*0.1


δγ    = 0.0*(exp(im*0))
γ     = γ + δγ

δμx   = -0.00*(U/8)
μx    = μx + δμx

Ω     = (μ0 .- U*U/(4.0*γ) .+ (γ*μx*μx)^(1.0/3.0)*ω)


rcParams["markers.fillstyle"] = "none"
hλ          = figure(num=1,figsize=[8.,6.]);
ax1         = gca()
pΛ          = ax1.plot(imag.(Ω),real.(Ω),linestyle="none",marker="o",markersize=8)
Ωi_max      = maximum(abs.(imag.(Ω)))
Ωr_min      = minimum(real.(Ω))
Ωr_max      = maximum(real.(Ω))
ax1.set_xlim((0.0,1.2*Ωi_max))
ax1.set_ylim((1.1*Ωr_min,0.5))
ax1.set_xlabel(L"\mathfrak{Im}(ω)",fontsize=lafs)
ax1.set_ylabel(L"\mathfrak{R}(ω)",fontsize=lafs)
# ax1.grid()
rcParams["markers.fillstyle"] = "full"

# Include the function files
include("sem_main.jl")

rng   = MersenneTwister(1235)

xg    = QT*(vimult.*Geom.xm1[:])
vt    = VT

#v     = randn(vt,ndof)*1.0;
v     = 1.0*exp.(-(xg.-5.75).^2)
v     = v .+ conj.(v)

ifplot      = true 
verbose     = true
nsteps      = 10000000
ifsave      = false
plotstep    = 20000
verbosestep = 20000
histstep    = 1000
nhist       = Int(nsteps/histstep)

Hist        = zeros(vt,nhist)
Time        = zeros(Float64,nhist)

#if (ifadjoint)
#  Ω = conj.(Ω)
#end  
if (ifplot)
  hv  = figure(num=2,figsize=[8.,6.]);
  ax2 = gca()
  ax2.set_xlabel(L"x",fontsize=lafs)
  ax2.set_ylabel(L"A",fontsize=lafs)
end

cm    = get_cmap("tab10");
rgba0 = cm(0) 
rgba1 = cm(1) 
rgba2 = cm(2) 

dt    = prec(0.001)

bc    = zeros(vt,1,ndof);
NGL(x)= NLGinzburgLandau(OPg,ones(vt,ndof),x,δ5,vt(0),vt(0),true,false)  

t = zro*dt        # Time
i = 0             # Istep

maxouter_it = 150
major_it    = 1


ifconv = false
println("Press x to stop. Any other key to continue...")
xin = readline()
if xin != "x"
  ifconv = false
else
  ifconv = true
end  

# Start iterations
println("Starting Iterations")

while (~ifconv)
  global v
  global t, i
  global plr,pli
  global OPg
  global hλ, hv, ax2

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

  if (mod(i,histstep) == 0)
    j = Int(i/histstep)
    Hist[j] = v[50]
    Time[j] = t
  end  

# Plot the field  
  if (ifplot && mod(i,plotstep)==0)
    if (i>plotstep) 
      for lo in ax2.get_lines()
        lo.remove()
      end  
    end  
   
    pv1 = ax2.plot(xg,real.(v),linestyle="-",color=rgba0)
    pv2 = ax2.plot(xg,imag.(v),linestyle="--",color=rgba0)
    pv2 = ax2.plot(xg,abs.(v) ,linestyle="-",color=rgba1,linewidth=2)

    vmax = 1.2*maximum(abs.(v))
    vmin = -vmax
    dv   = abs(vmax-vmin)
    ax2.set_ylim((vmin,vmax))
#    ax2.set_ylim((-dv,dv))
    hv.show()    
   
  end 
  
  if i==nsteps
    break
  end  

end       # while ... 

if !ifconv
  h3  = figure(num=3,figsize=[8.,6.]);
  ax3 = gca()
  ax3.plot(Time,real.(Hist))
  ax3.set_xlabel(L"t",fontsize=lafs)
  ax3.set_ylabel(L"A",fontsize=lafs)
  
  
#  # Frequency
#  peak_ind    = argmaxima(real.(Hist))
#  peak_times  = Time[peak_ind]
#  delta_times = diff(peak_times)
#  afreq       = 2.0*π./delta_times
#  npeaks      = length(afreq)
#  if npeaks > 9
#    ωend      = afreq[npeaks-9:npeaks]
#  else
#    ωend      = afreq
#    println("Not enough cycles: $npeaks")
#  end
#  
#  ωmean       = mean(ωend)
#  @printf("δμx : %.5f ; Ω: %.4e \n", δμx, ωmean)
#  @printf("δγ_r: %.5f ; Ω: %.4e \n", real(δγ), ωmean)
#
#  # Saved values
#  δμx_v     = zeros(Float64,7)
#  Ω_v1      = zeros(Float64,7)
#
#  δμx_v[1]  = -0.00313
#  Ω_v1[1]   =  1.0139
# 
#  δμx_v[2]  = -0.00625
#  Ω_v1[2]   =  1.0272
#
#  δμx_v[3]  = -0.0125
#  Ω_v1[3]   =  1.0519
#
#  δμx_v[4]  = -0.01875
#  Ω_v1[4]   =  1.0753
#
#  δμx_v[5]  = -0.025
#  Ω_v1[5]   =  1.0977
#
#  δμx_v[6]  = -0.03125
#  Ω_v1[6]   =  1.1194
#
#  δμx_v[7]  = -0.0375
#  Ω_v1[7]   =  1.1403
#
#  δγ_v      = zeros(Float64,10)
#  Ω_v2      = zeros(Float64,10)
#
#  δγ_v[1]   = -0.0250
#  Ω_v2[1]   =  1.0037
#
#  δγ_v[2]   = -0.0500
#  Ω_v2[2]   =  1.0081
 
end


if (ifsave)
  fname = "GL_BaseState5.jld2"
  save(fname,"xg",xg,"basestate",v,"U",U,"γ",γ,"μx",μx,"μ0",μ0,"δ5",δ5,"δγ",δγ,"δμx",δμx);
  println(fname*" saved.")
end 

println("Done.")












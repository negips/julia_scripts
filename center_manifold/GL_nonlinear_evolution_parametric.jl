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

ω  = [ω1, ω2, ω3, ω4, ω5, ω6, ω7, ω8, ω9, ω10, ω11, ω12, ω13, ω14, ω15]
#Ω  = im*(U*U/8.0 .- U*U/(4.0*γ) .+ γ^(1.0/3.0)*(U^(4.0/3.0))/(160.0^(2.0/3.0))*ω)

# Ginzburg Landau Parameters
Ω0    = im*1.0
U     = 1.0
ϕ     = -π/4.0
R     = 1.0

γ_0   = R*exp(im*ϕ)
γ     = copy(γ_0)

μx_0  = U/8.0 
μx    = copy(μx_0)

μ0    = Ω0 + (U^2)/(4.0*γ) - ((γ*μx*μx)^(1.0/3.0))*ω1
δ5    = (-1.0 + 1.0im)*0.1

Ω     = (μ0 .- U*U/(4.0*γ) .+ (γ*μx*μx)^(1.0/3.0)*ω)

nγ    = 1
if nγ == 1
  δγ_v  = zeros(Float64,1)
else
  δγ_v  = LinRange(0.0,0.25,nγ)*abs(γ_0)
end  


nμx   = 10
if nμx == 1
  δμx_v = zeros(Float64,1) 
else
  δμx_v = LinRange(0.0,0.25,nμx)*abs(μx_0)
end  

ifplot      = true
verbose     = false
nsteps      = 200000
ifsave      = false
plotstep    = 2000
verbosestep = 2000
histstep    = 10
nhist       = Int(nsteps/histstep)

if (ifplot)
  hv  = figure(num=2,figsize=[8.,6.]);
  ax2 = gca()
  ax2.set_xlabel(L"x",fontsize=lafs)
  ax2.set_ylabel(L"A",fontsize=lafs)
end

ω_matrix    = zeros(Float64,nγ,nμx)

for iγ in 1:nγ
for iμ in 1:nμx

  global hv
  global ax2
  global ω_matrix
  global γ, μx
  
  δγ  = -δγ_v[iγ]
  δμx = -δμx_v[iμ]

  @printf("\nδμx : %.5f ; δγ_r: %.5f \n\n", δμx, real(δγ))

  γ   = γ_0  + δγ
  μx  = μx_0 + δμx

  if abs(δγ)< 100.0*eps(1.0) && abs(δμx) < 100.0*eps(1.0)
    ω_matrix[1,1] = imag(Ω0)
    continue
  end

  # Include the function files
  include("sem_main.jl")

  rng   = MersenneTwister(1235)
  
  xg    = QT*(vimult.*Geom.xm1[:])
  vt    = VT
  
  v     = randn(vt,ndof)*1.0;
  
  
  Hist  = zeros(vt,nhist)
  Time  = zeros(Float64,nhist)
  
  cm    = get_cmap("tab10");
  rgba0 = cm(0) 
  rgba1 = cm(1) 
  rgba2 = cm(2) 
  
  dt    = prec(0.001)
  
  bc    = zeros(vt,1,ndof);
  NGL(x)= NLGinzburgLandau(OPg,ones(vt,ndof),x,δ5,vt(0),vt(0),true,false)  

  t = zro*dt        # Time


  # Start iterations
  println("Starting Iterations")

  for i in 1:nsteps
  
    t = t + dt;
  
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
      hv.show()    
    end 
  
  end       # i in 1:nsteps 

  
  
  # Frequency
  peak_ind    = argmaxima(real.(Hist))
  peak_times  = Time[peak_ind]
  delta_times = diff(peak_times)
  afreq       = 2.0*π./delta_times
  npeaks      = length(afreq)
  if npeaks > 9
    ωend      = afreq[npeaks-9:npeaks]
  else
    ωend      = afreq
    println("Not enough cycles: $npeaks")
  end
  
  ωmean       = mean(ωend)
  @printf("\nδμx : %.5f ; δγ_r: %.5f ; Ω: %.4e \n\n", δμx, real(δγ), ωmean)
  # @printf("δγ_r: %.5f ; Ω: %.4e \n", real(δγ), ωmean)
  ω_matrix[iγ,iμ] = ωmean

end   # iμ  
end   # iγ  

println("Done.")












println("Non-Linear evolution for Ginzburg Landau equations")

#include("Module_SEM1D/SEM1D.jl")
#using .SEM1D

#include("Module_StepperArnoldi/StepperArnoldi.jl")
##using .StepperArnoldi

#include("Module_CenterManifold/CenterManifold.jl")
##using .CenterManifold

using Peaks
using Statistics
using Random
using JLD2

include("$JULIACOMMON/RK4.jl")
include("NLGinzburgLandau.jl")
include("OP_RK4.jl")


lafs        = 16
lgfs        = 12
mksz        = 6

include("GL_Setup.jl")
#-------------------------------------------------- 

screen            = 2

if screen == 1
  # hp spectre
  lafs      = 16        # Label font size
  lgfs      = 12        # Legend font size
  figsz1    = [8.0, 6.0]
  figsz2    = [17.0, 6.0]
elseif screen == 2
  # workstation
  lafs      = 20        # Label font size
  lgfs      = 16        # Legend font size
  figsz1    = [24.0, 20.0]
  figsz2    = [34.0, 20.0]
end  

close("all")

ifplot      = true
histplot    = true
verbose     = true
nsteps      = 30000000
ifsave      = false
plotstep    = 20000
verbosestep = 20000
histstep    = 1000
nhist       = Int(nsteps/histstep)
hist_x      = 10.0                        # Location of history point
hist_i      = argmin(abs.(xg .- hist_x))  # Index of history point
nfreq       = 1                           # No. of external frequencies
dt          = 0.0001

θA          = [0.1] #; 0.1; 0.25; 0.5; 1.0]
nθ          = length(θA)


cm          = get_cmap("tab10");
rgba0       = cm(0) 
rgba1       = cm(1) 
rgba2       = cm(2) 


vt          = Complex{Inp.Dtype}
zro         = vt(0)

TLast       = 2500.0
Hist        = zeros(vt,nhist,nθ)
Time        = zeros(Float64,nhist)
Peak_Amp    = zeros(Float64,nθ)
ω_nonlinear = zeros(Float64,nθ)


# Work Arrays
vwork       = zeros(vt,ndof,5)
θwork       = zeros(vt,nfreq,5)

if (ifplot)
  hv  = figure(num=2,figsize=figsz1);
  ax2 = gca()
  ax2.set_xlabel(L"x",fontsize=lafs)
  ax2.set_ylabel(L"A",fontsize=lafs)
end

if (histplot)
  h3  = figure(num=3,figsize=figsz2);
  ax3 = gca()
end  

NGL(x)= NLGinzburgLandau(OPg,ones(vt,ndof),x,δ[5],zro,zro,Inp.lbc,Inp.rbc)  

# Forcing Shape
x0    = ForcingLocation()
ψ     = zeros(ComplexF64,ndof)
SetForcingShape!(ψ,Bg,xg,x0,1.0)
F     = zeros(ComplexF64,ndof,nfreq)
copyto!(F,ψ)
Ω     = [1.3im]
FNGL(x,y) = ForcedNLGinzburgLandau(OPg,ones(vt,ndof),x,F,y,δ[5],Ω,zro,zro,Inp.lbc,Inp.rbc)

println("Press x to stop. Any other key to continue...")
xin = readline()
if xin == "x"
  nsteps = 0
end  
# Start iterations
println("Starting Iterations")

for ik in 1:nθ
  global vwork,θwork
  global hv, ax2
  global h3, ax3

  θAmp  = θA[ik]
  t     = Inp.Dtype(0)    # Time
  
  rng   = Xoshiro(1235)
  rnd   = rand(rng,vt,ndof)
  z     = zeros(vt,6)
  z[1]  = 1.0e-5*vt(1)
  z[2]  = z[1]'
  z[5]  = θAmp*vt(1)
  z[6]  = z[5]'
  fld   = Get_AsymptoticField(z,Vext,Y_O2,Y_O3)

  v     = zeros(vt,ndof)
  # v     = fld[1:ndof]

  #v    .= v .+ 1.0e-5*rnd
  #v    .= v .+ 0.1*exp.(-(xg.-5.75).^2)
  #v    .= v .+ conj.(v)
  θ     = zeros(vt,nfreq)
  θ[1]  = θAmp*(1.0 + 0.0im)

  # Testing temporary forcing amplitude change
  θtmp  = 1.0
  λtmp  = -0.01

  for i in 1:nsteps
  
    t = t + dt;
  
    # Apply BC
    SEM1D.SEM_SetBC!(v,Inp.lbc,Inp.rbc)
    # Non-linear Evolution
    # OP_RK4!(NGL,v,dt)

    # Testing temporary forcing amplitude change
    θ = θ .+ θtmp*exp(λtmp*t)  

    # Forced Non-linear Evolution
    # OP2_RK4!(FNGL,v,θ,dt)
    OP2_RK4!(FNGL,v,θ,dt,vwork,θwork)

    # Testing temporary forcing amplitude change
    θ = θ .- θtmp*exp(λtmp*t)  
   
    # Print something  
    if verbose && mod(i,verbosestep)==0
      println("Istep=$i, Time=$t")
    end
  
    if (mod(i,histstep) == 0)
      j = Int(i/histstep)
      Hist[j,ik]  = v[hist_i]
      Time[j]     = t
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
      # ax2.set_ylim((-dv,dv))
      hv.show()    
    end 

    # History plot
    if histplot && mod(i,plotstep)==0
      if (i>plotstep) 
        for lo in ax3.get_lines()
          lo.remove()
        end  
      end  
      ax3.plot(Time[1:j],real.(Hist[1:j,ik]),color=cm(ik-1))
    end

  end       # i in 1:nsteps

  if nsteps>0 && histplot
    ax3.plot(Time,real.(Hist[:,ik]))
    ax3.set_xlabel(L"t",fontsize=lafs)
    ax3.set_ylabel(L"A",fontsize=lafs)
 
    # ax3.set_xlim([2900.0,3000.0])

    linds           = Time .> TLast
    time2           = Time[linds]
    hist2           = Hist[linds,:]
    pkind           = argmaxima(real.(hist2[:,ik]))
    pktimes         = time2[pkind]
    mamp            = real.(hist2[pkind,ik])
    Peak_Amp[ik]    = mean(mamp) 
    delta_times     = diff(pktimes)
    afreq           = 2.0*π./delta_times
    ω_nonlinear[ik] = mean(afreq)

    @printf("|θ|: %.2f ; Amax: %.5f ; Ω: %.4e\n", θAmp,Peak_Amp[ik], ω_nonlinear[ik])
  end       # if nsteps>0 && histplot

end         # ik in 1:nθ

ax3.set_xlim([2900.0,3000.0])


# Plot Peaks
if nsteps>0 && histplot
  h4          = figure(num=4,figsize=[8.,6.]);
  ax4         = gca()
  last_inds   = Time .> TLast
  Time2       = Time[last_inds]
  Hist2       = Hist[last_inds,:]
  Peak_Amp    = zeros(Float64,nθ)
  for ik in 1:nθ
    peak_ind        = argmaxima(real.(Hist2[:,ik]))
    peak_times      = Time2[peak_ind]
    maxamp          = real.(Hist2[peak_ind,ik])
    Peak_Amp[ik]    = mean(maxamp) 
    delta_times     = diff(peak_times)
    afreq           = 2.0*π./delta_times
    ω_nonlinear[ik] = mean(afreq)
  end
  ax4.plot(θA,Peak_Amp,linestyle="none",marker="o",markersize=mksz)
end  



if (ifsave)
  # fname = "GL_BaseState5.jld2"
  # save(fname,"xg",xg,"basestate",v,"U",U,"γ",γ,"μx",μx,"μ0",μ0,"δ5",δ5,"δγ",δγ,"δμx",δμx);
  # println(fname*" saved.")
end 

println("Done.")












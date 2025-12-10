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


include("GL_Setup.jl")
#-------------------------------------------------- 

close("all")

ifplot      = true
histplot    = true
moveaxis    = true
verbose     = true
if ifresonant
  nsteps    = 30000000
else
  nsteps    = 30000000
end  

ifsave      = true
plotstep    = 10000
verbosestep = 10000
histstep    = 1000
nhist       = Int(nsteps/histstep)
hist_x      = ForcingLocation()           # Location of history point
hist_i      = argmin(abs.(xg .- hist_x))  # Index of history point
nfreq       = 1                           # No. of external frequencies
dt          = 0.0001
Tend        = dt*nsteps

#θA          = [0.1; 0.25; 0.5; 0.75; 1.0]
θA          = [0.1; 0.2; 0.3; 0.4; 0.5]
#θA          = [0.5]
nθ          = length(θA)
ncycles     = ones(Int64,nθ)
if !ifresonant
  ncycles[1]  = 3
  ncycles[2]  = 2
end

cm          = get_cmap("tab10");
rgba0       = cm(0) 
rgba1       = cm(1) 
rgba2       = cm(2) 


vt          = Complex{Inp.Dtype}
zro         = vt(0)

TLast       = Tend - 500.0
Hist        = zeros(vt,nhist,nθ)
Time        = zeros(Float64,nhist)
Peak_Amp    = zeros(Float64,nθ)
ω_nonlinear = zeros(Float64,nθ)
vlast       = zeros(vt,ndof,nθ)


# Work Arrays
vwork       = zeros(vt,ndof,5)
θwork       = zeros(vt,nfreq,5)

if (ifplot)
  hv  = figure(num=2,figsize=Grh.figsz1);
  ax2 = gca()
  ax2.set_xlabel(L"x",fontsize=Grh.lafs)
  ax2.set_ylabel(L"A",fontsize=Grh.lafs)
end

if (histplot)
  h3  = figure(num=3,figsize=Grh.figsz2);
  ax3 = gca()
end  

NGL(x)= NLGinzburgLandau(OPg,ones(vt,ndof),x,δ[5],zro,zro,Inp.lbc,Inp.rbc)  

# Forcing Shape
x0,κ  = ForcingParams()
ψ     = zeros(ComplexF64,ndof)
SetForcingShape!(ψ,Bg,xg,x0,1.0,κ)
F     = zeros(ComplexF64,ndof,nfreq)
copyto!(F,ψ)
Ωf    = [λh[1]]
FNGL(x,y) = ForcedNLGinzburgLandau(OPg,ones(vt,ndof),x,F,y,δ[5],Ωf,zro,zro,Inp.lbc,Inp.rbc)
Vdamp     = zeros(vt,ndof,1)
Vdamp[:,1]= V[ind1,1]
Wdamp     = zeros(vt,ndof,1)
Wdamp[:,1]= W[ind1,1]
χdamp     = [vt(0.00/dt)]
DFGL(x,y) = DampedFGL(FNGL,Bg,x,y,Vdamp,Wdamp,χdamp)

# For θ evolution
ΩM        = zeros(vt,1,1)
ΩM[1,1]   = Ωf[1]
FΩ(x)     = StuartLandau1(ΩM,x)

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
  global vlast

  θAmp  = θA[ik]
  
  rng   = Xoshiro(1235)
  rnd   = rand(rng,vt,ndof)
  z     = zeros(vt,6)
  z[1]  = 0.0e-1*vt(1)
  z[2]  = z[1]'
  z[5]  = θAmp*vt(1)
  z[6]  = z[5]'
  fld   = Get_AsymptoticField(z,Vext,Y_O2,Y_O3)

  v     = zeros(vt,ndof)
  #v     = fld[1:ndof]

  #v    .= v .+ 1.0e-5*rnd
  #v    .= v .+ 0.1*exp.(-(xg.-5.75).^2)
  #v    .= v .+ conj.(v)
  global θ     = zeros(vt,nfreq)
  θ[1]  = θAmp*(1.0 + 0.0im)

  # Testing temporary forcing amplitude change
  θtmp  = vt(1.00)
  λtmp  = -0.025

  cycles = ncycles[ik]
  println("$cycles cycles for ik=$ik")
  for ic = 1:cycles

    t     = Inp.Dtype(0)    # Time

    for i in 1:nsteps
    
      t = t + dt;
    
      # Apply BC
      SEM1D.SEM_SetBC!(v,Inp.lbc,Inp.rbc)
      # Non-linear Evolution
      # OP_RK4!(NGL,v,dt)

      # Testing temporary forcing amplitude change
      if ic==1
        fac  = (1.0 - exp(λtmp*t))
      else
        fac  = 1.0
      end
      θ    = θ*fac

      # Forced Non-linear Evolution
      # OP2_RK4!(FNGL,v,θ,dt)
      # OP2_RK4!(FNGL,v,θ,dt,vwork,θwork)
      if (ik<2 && ic==1)
        OP2_RK4!(DFGL,v,θ,dt,vwork,θwork)
      else
        OP2_RK4!(FNGL,v,θ,dt,vwork,θwork)
      end  


      # Testing temporary forcing amplitude change
      #OP_RK4!(FΩ,dθ,dt)
      #θ = θ .- dθ
      θ    = θ/fac

      χd   = abs.(Wdamp'*(Bg.*v))[1]

      # Print something  
      if verbose && mod(i,verbosestep)==0
        println("ik=$ik/$nθ, ic=$ic/$cycles, Istep=$i, Time=$t, θ=$(abs.(θ)[1]), χd=$(χd)")
      end
   
      # Save History
      if (mod(i,histstep) == 0)
        j = Int(i/histstep)
        Hist[j,ik]  = v[hist_i]
        Time[j]     = t
      end  
    
      # Plot the field  
      if (ifplot && mod(i,plotstep)==0)
        #if (i>plotstep) 
          for lo in ax2.get_lines()
            lo.remove()
          end  
        #end  
       
        pv1 = ax2.plot(xg,real.(v),linestyle="-",color=rgba0)
        pv2 = ax2.plot(xg,imag.(v),linestyle="--",color=rgba0)
        pv2 = ax2.plot(xg,abs.(v) ,linestyle="-",color=rgba1,linewidth=2)
    
        vmax = 1.2*maximum(abs.(v))
        vmin = -vmax
        dv   = abs(vmax-vmin)
        ax2.set_ylim((vmin,vmax))
        # ax2.set_ylim((-dv,dv))
        hv.show()    

        # History plot
        if histplot
          #if (i>plotstep) 
            for lo in ax3.get_lines()
              lo.remove()
            end  
          #end  
          ax3.plot(Time[1:j],real.(Hist[1:j,ik]),color=cm(ik-1))
        end

        if (moveaxis)
          tmax = Time[j]
          tmin = max(0.0,tmax-300.0)
          ax3.set_xlim([tmin,tmax])
        end  
      end   # ifplot 
    end     # i in 1:nsteps
  end       # ic in 1:cycles    

  vlast[:,ik] = copy(v)
  if nsteps>0 && histplot

    # Remove previous plots
    for lo in ax3.get_lines()
      lo.remove()
    end  
  
    # Plot entire history
    ax3.plot(Time,real.(Hist[:,ik]))
    ax3.set_xlabel(L"t",fontsize=Grh.lafs)
    ax3.set_ylabel(L"A",fontsize=Grh.lafs)
 
    # ax3.set_xlim([2900.0,3000.0])
    linds           = Time .> TLast
    time2           = Time[linds]
    hist2           = Hist[linds,ik]
    pkind           = argmaxima(real.(hist2))
    pktimes         = time2[pkind]
    mamp            = real.(hist2[pkind])
    Peak_Amp[ik]    = mean(mamp) 
    delta_times     = diff(pktimes)
    afreq           = 2.0*π./delta_times
    ω_nonlinear[ik] = mean(afreq)

    @printf("|θ|: %.2f ; Amax: %.5f ; Ω: %.4e\n", θAmp,Peak_Amp[ik], ω_nonlinear[ik])
  end       # if nsteps>0 && histplot
end         # ik in 1:nθ

# ax3.set_xlim([2900.0,3000.0])


# Plot Peaks
if nsteps>0 && histplot
  h4          = figure(num=4,figsize=Grh.figsz1);
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

    @printf("|θ|: %.2f ; Amax: %.5f ; Ω: %.4e\n", θA[ik],Peak_Amp[ik], ω_nonlinear[ik])
  end
  ax4.plot(θA,Peak_Amp,linestyle="none",marker="o",markersize=Grh.mksz)
end  


if (ifsave && nsteps>0)
  if ifresonant
    fname = "GL_resonant_Parametric3.jld2"
  else
    fname = "GL_nonresonant_Parametric3.jld2"
  end  
  save(fname,"xg",xg,"vlast",vlast,"Vext",Vext,"Y_O2",Y_O2,"Y_O3",Y_O3,"Khat",Khat,"G_O2",G_O2,"G_O3",G_O3,"δ",δ,"Time",Time,"θA",θA,"Peak_Amp",Peak_Amp,"Hist",Hist,"ω_nonlinear",ω_nonlinear);
  println(fname*" saved.")
end 

println("Done.")












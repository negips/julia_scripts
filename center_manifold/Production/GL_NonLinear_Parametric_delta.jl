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


include("GL_Setup2.jl")
#-------------------------------------------------- 

close("all")

ifplot      = true
histplot    = true
moveaxis    = true
verbose     = true
nsteps      = 30000000

ifsave      = true
plotstep    = 10000
verbosestep = 10000
histstep    = 1000
nhist       = Int(nsteps/histstep)
hist_x,tmp  = ForcingParams()             # Location of history point
hist_i      = argmin(abs.(xg .- hist_x))  # Index of history point
nfreq       = 1                           # No. of external frequencies
dt          = 0.0001
Tend        = dt*nsteps

δ4          = -Vector(1:10)*0.02
nδ4         = length(δ4)
ncycles     = ones(Int64,nδ4)
for i in 1:3
  ncycles[i]  = 2
end
for i in 4:4
  ncycles[i]  = 1
end

cm          = get_cmap("tab10");
rgba0       = cm(0) 
rgba1       = cm(1) 
rgba2       = cm(2) 


vt          = Complex{Inp.Dtype}
zro         = vt(0)

TLast       = Tend - 500.0
Hist        = zeros(vt,nhist,nδ4)
Time        = zeros(Float64,nhist)
Peak_Amp    = zeros(Float64,nδ4)
ω_nonlinear = zeros(Float64,nδ4)
vlast       = zeros(vt,ndof,nδ4)

# Work Arrays
vwork       = zeros(vt,ndof,5)

if (ifplot)
  hv  = figure(num=2,figsize=Grh.figsz3);
  ax2 = gca()
  ax2.set_xlabel(L"x",fontsize=Grh.lafs)
  ax2.set_ylabel(L"A",fontsize=Grh.lafs)
end

if (histplot)
  h3  = figure(num=3,figsize=Grh.figsz2);
  ax3 = gca()
end  

# NGL(x)= NLGinzburgLandau(OPg,Bg,x,δ[5],zro,zro,Inp.lbc,Inp.rbc)  


println("Press x to stop. Any other key to continue...")
xin = readline()
if xin == "x"
  nsteps = 0
end  
# Start iterations
println("Starting Iterations")

for ik in 1:nδ4
  global vwork
  global hv, ax2
  global h3, ax3
  global vlast
  global δ, δc

 
  rng   = Xoshiro(1235)
  v     = 1.0e-4*ones(vt,ndof)

  # GL Parameters
  δ4_diff   = δ4[ik]
  δ         = Set_GL_CriticalParams()
  δ[4]      = δ[4] + δ4_diff
  δc        = conj.(δ)

  # Reset Operators
  include("GL_OP_Setup.jl")
  NGL(x) = NLGinzburgLandau(OPg,Bg,x,δ[5],zro,zro,Inp.lbc,Inp.rbc)  

  cycles = ncycles[ik]
  println("$cycles cycles for ik=$ik")
  for ic = 1:cycles

    t     = Inp.Dtype(0)    # Time

    for i in 1:nsteps
    
      t = t + dt;
    
      # Apply BC
      SEM1D.SEM_SetBC!(v,Inp.lbc,Inp.rbc)

      # Non-linear Evolution
      OP_BiRK4!(NGL,Bgi,v,dt,vwork)

      # Print something  
      if verbose && mod(i,verbosestep)==0
        println("ik=$ik/$nδ4, ic=$ic/$cycles, Istep=$i, Time=$t")
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
        pv2 = ax2.plot(xg,abs.(v) ,linestyle="-",color=rgba1,linewidth=3)
    
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

    @printf("δ4: %.2f ; Amax: %.5f ; Ω: %.4e\n", δ4_diff,Peak_Amp[ik], ω_nonlinear[ik])
  end       # if nsteps>0 && histplot
end         # ik in 1:nδ4

# ax3.set_xlim([2900.0,3000.0])


# Plot Peaks
if nsteps>0 && histplot
  h4          = figure(num=4,figsize=Grh.figsz1);
  ax4         = gca()
  last_inds   = Time .> TLast
  Time2       = Time[last_inds]
  Hist2       = Hist[last_inds,:]
  Peak_Amp    = zeros(Float64,nδ4)
  for ik in 1:nδ4
    peak_ind        = argmaxima(real.(Hist2[:,ik]))
    peak_times      = Time2[peak_ind]
    maxamp          = real.(Hist2[peak_ind,ik])
    Peak_Amp[ik]    = mean(maxamp) 
    delta_times     = diff(peak_times)
    afreq           = 2.0*π./delta_times
    ω_nonlinear[ik] = mean(afreq)

    @printf("δ4: %.2f ; Amax: %.5f ; Ω: %.4e\n", δ4[ik], Peak_Amp[ik], ω_nonlinear[ik])
  end
  ax4.plot(-δ4,ω_nonlinear,linestyle="none",marker="o",markersize=Grh.mksz)
end  


if (ifsave && nsteps>0)
  fname = "GL_diffusion_Parametric2.jld2"
  δ         = Set_GL_CriticalParams()

  save(fname,"xg",xg,"vlast",vlast,"δ",δ,"Time",Time,"δ4",δ4,"Peak_Amp",Peak_Amp,"Hist",Hist,"ω_nonlinear",ω_nonlinear);
  println(fname*" saved.")
end 

println("Done.")












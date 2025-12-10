println("Non-Linear evolution for the Stuart Landau equations")

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

# include("GL_Setup.jl")
#-------------------------------------------------- 

close("all")

cm          = get_cmap("tab10");
rgba0       = cm(0) 
rgba1       = cm(1) 
rgba2       = cm(2) 

ifplot      = true
histplot    = true
moveaxis    = true
plotfield   = true
verbose     = true
if ifresonant
  nsteps    = 3000000
else
  nsteps    = 3000000
end
ifsave      = false
plotstep    = 20000
verbosestep = 10000
histstep    = 100
nhist       = Int(nsteps/histstep)
xhist       = true
vt          = Complex{Inp.Dtype}
zro         = vt(0)

dt          = 0.001
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

Hist_Mode   = zeros(vt,nhist,m,nθ)
Time        = zeros(Float64,nhist)
Peak_Amp    = zeros(Float64,nθ)
ω_nonlinear = zeros(Float64,nθ)
Mode_Ind    = [1]                         # Which mode to plot 

figsz       = [12.0, 5.0]

h3          = figure(num=3,figsize=figsz);
ax3         = gca()
ax3.set_xlabel(L"time",fontsize=lafs)
ax3.set_ylabel(L"Z_{i}",fontsize=lafs)

TLast       = Tend - 500.0

if xhist
  hist_x    = ForcingLocation()           # Location of history point
  hist_i    = argmin(abs.(xg .- hist_x))  # Index of history point

  Histx     = zeros(vt,nhist,nθ)

  h4        = figure(num=4,figsize=figsz);
  ax4       = gca()
  ax4.set_xlabel(L"time",fontsize=lafs)
  ax4.set_ylabel(L"A_{x}",fontsize=lafs)
end  

if plotfield
  h5        = figure(num=5,figsize=figsz);
  ax5       = gca()
  ax5.set_xlabel(L"x",fontsize=lafs)
  ax5.set_ylabel(L"A",fontsize=lafs)
end

# Stuart Landau
G1          = Khat
G2          = G_O2
G3          = G_O3
SL(x)       = StuartLandau3(G1,G2,G3,x)  

println("Press x to stop. Any other key to continue...")
xin = readline()
if xin == "x"
  nsteps = 0
end  


println("Starting Iterations")

for ik in 1:nθ

  # Initialize z
  θAmp        = θA[ik]

  z           = zeros(vt,m)
  rng         = Xoshiro(1235)
  # Mode initial values
  z[1]        = 1.0e-4*rand(rng,vt)
  z[2]        = z[1]'
  # Parameter Perturbations
  z[n+1:n+p]  = zeros(vt,p)
  # Harmonic Forcing Amplitude
  z[n+p+1]    = θAmp*(1.0 + 0.0im)
  z[n+p+2]    = z[n+p+1]'
  
  # Work Arrays
  zwork       = zeros(vt,m,5)
  

  # Testing temporary forcing amplitude change
  θtmp  = vt(1.00)
  λtmp  = -0.01

  cycles = ncycles[ik]
  for ic in 1:cycles

    # Start iterations
    t   = Inp.Dtype(0)    # Time
    
    for i in 1:nsteps
    
      t = t + dt;
 
      # Testing temporary forcing amplitude change
      # if ic == 1
      #   fac  = (1.0 - exp(λtmp*t))
      # else
      #   fac  = 1.0
      # end
      # z[n+p+1]    = z[n+p+1]*fac
      # z[n+p+2]    = z[n+p+2]*fac

      # Stuart Landau Evolution
      OP_RK4!(SL,z,dt,zwork)

      # z[n+p+1]    = z[n+p+1]/fac
      # z[n+p+2]    = z[n+p+2]/fac

      # Set conjugation correctly
      z[2] = z[1]'
    
      znorm = sqrt(abs(z'*z))
      # Print something  
      if verbose && mod(i,verbosestep)==0
        println("ik=$ik/$nθ, ic=$ic/$cycles, Istep=$i, Time=$t, |z|=$znorm")
      end
    
      if (mod(i,histstep) == 0)
        j = Int(i/histstep)
        Hist_Mode[j,:,ik] = copy(z)
        Time[j]           = t
    
        # Get field value at point x = hist_x,
        # corresponding to array index hist_i
        if (xhist)
          Histx[j,ik] = Get_AsymptoticFieldx(hist_i,z,Vext,Y_O2,Y_O3)
        end  
      end

      if ifplot && mod(i,plotstep)==0

        # Remove previous plots
        for lo in ax3.get_lines()
          lo.remove()
        end  
        ax3.plot(Time[1:j],real.(Hist_Mode[1:j,Mode_Ind,ik]),color=cm(ik-1))

        # Remove previous plots
        for lo in ax4.get_lines()
          lo.remove()
        end
        ax4.plot(Time[1:j],real.(Histx[1:j,ik]),color=cm(ik-1))

        if (moveaxis)
          tmax = Time[j]
          tmin = max(0.0,tmax-500.0)
          ax3.set_xlim([tmin,tmax])
          ax4.set_xlim([tmin,tmax])
        end  

        if (plotfield)
          # Remove previous plots
          for lo in ax5.get_lines()
            lo.remove()
          end
          fld12 = CenterManifold.GetAsymptoticField3(z,Vext,Y_O2,Y_O3)
          fld1  = fld12[1:ndof]
          ax5.plot(xg,real.(fld1),color=cm(0),linestyle="-",linewidth=1)
          ax5.plot(xg,imag.(fld1),color=cm(0),linestyle="--",linewidth=1)
          ax5.plot(xg,abs.(fld1),color=cm(0),linestyle="-",linewidth=3)
        end  
      end   # ifplot && mod(i,plotstep)==0 
    end     # i in 1:nsteps
  end       # ic in 1:cycles

  if histplot && nsteps>0
    # Remove previous plots
    for lo in ax3.get_lines()
      lo.remove()
    end  
    ax3.plot(Time,real.(Hist_Mode[:,Mode_Ind,ik]))
    # ax3.set_xlabel(L"t",fontsize=lafs)
    # ax3.set_ylabel(L"A",fontsize=lafs)

    # Remove previous plots
    for lo in ax4.get_lines()
      lo.remove()
    end  
    ax4.plot(Time,real.(Histx[:,ik]),color=cm(ik-1))

    linds           = Time .> TLast
    time2           = Time[linds]
    hist2           = Histx[linds,ik]
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

#ax4.set_xlim([4500.0,5000.0])
#ax4.set_ylim([-0.4,0.4])
#ax4.legend(fontsize=lgfs,ncols=nθ)


# Plot Peaks
if nsteps>0 && histplot
  h6          = figure(num=6,figsize=[8.,6.]);
  ax6         = gca()
  last_inds   = Time .> TLast
  Time2       = Time[last_inds]
  Hist2       = Histx[last_inds,:]
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
  ax6.plot(θA,Peak_Amp,linestyle="none",marker="o",markersize=mksz)
  ax6.set_xlabel(L"|θ|",fontsize=lafs)
  ax6.set_ylabel(L"A_{x}^{max}",fontsize=lafs)
end  


if (ifsave && nsteps>0)
  if ifresonant
    fname = "SL_resonant_Parametric2.jld2"
  else
    fname = "SL_nonresonant_Parametric2.jld2"
  end
  save(fname,"xg",xg,"Vext",Vext,"Y_O2",Y_O2,"Y_O3",Y_O3,"G1",G1,"G2",G2,"G3",G3,"δ",δ,"Time",Time,"θA",θA,"Peak_Amp",Peak_Amp,"Histx",Histx,"ω_nonlinear",ω_nonlinear);
  println(fname*" saved.")
end 


println("Done.")












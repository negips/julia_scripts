println("Non-Linear evolution for the Stuart Landau equations")

using Peaks
using Statistics
using Random
using JLD2

include("$JULIACOMMON/RK4.jl")
include("NLGinzburgLandau.jl")
include("OP_RK4.jl")


# lafs        = 16
# lgfs        = 12
# mksz        = 6

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
#θA          = [0.1; 0.2; 0.3; 0.4; 0.5]
#θA          = Vector(1:10)*0.1
θA          = [0.25]
nθ          = length(θA)
ncycles     = ones(Int64,nθ)
if !ifresonant
  ncycles[1]  = 4
  ncycles[2]  = 3
end

Hist_Mode   = zeros(vt,nhist,m,nθ)
Time        = zeros(Float64,nhist)
Peak_Amp    = zeros(Float64,nθ)
ω_nonlinear = zeros(Float64,nθ)
Mode_Ind    = [1; 3]                         # Which mode(s) to plot 

# figsz       = [12.0, 5.0]

h3          = figure(num=3,figsize=Grh.figsz3);
ax3         = gca()
ax3.set_xlabel(L"time",fontsize=Grh.lafs)
ax3.set_ylabel(L"Z_{i}",fontsize=Grh.lafs)

TLast       = Tend - 500.0

if xhist
  hist_x,tmp= ForcingParams()             # Location of history point
  hist_i    = argmin(abs.(xg .- hist_x))  # Index of history point

  Histx     = zeros(vt,nhist,nθ)

  h4        = figure(num=4,figsize=Grh.figsz3);
  ax4       = gca()
  ax4.set_xlabel(L"time",fontsize=Grh.lafs)
  ax4.set_ylabel(L"A_{x}",fontsize=Grh.lafs)
end  

if plotfield
  h5        = figure(num=5,figsize=Grh.figsz3);
  ax5       = gca()
  ax5.set_xlabel(L"x",fontsize=Grh.lafs)
  ax5.set_ylabel(L"A",fontsize=Grh.lafs)
end

# Stuart Landau
G1          = Khat
G2          = G2
G3          = G3
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
  for i in 1:nsys
    if mod(i-1,2) == 0
      z[i]        = 1.0e-4*rand(rng,vt)
    else
      z[i]        = z[i-1]'
    end
  end
  # System Perturbations
  for i in 1:nsys
    j = PertModesExt[i]
    if j != 0
      z[j] = -σext[i]
    end
  end  

  # Parameter Perturbations
  for i in nsys+npert+1:nsys+npert+p
    z[i]  = 0
  end
  # Harmonic Forcing Amplitude
  for i in nsys+npert+p+1:m
    j = i - (nsys+npert+p)
    if mod(j-1,2) == 0 
      z[i]    = θAmp*(1.0 + 0.0im)
    else  
      z[i]    = z[i-1]'
    end
  end
  # println(z)

  # Work Arrays
  zwork       = zeros(vt,m,5)

  # Testing temporary forcing amplitude change
  θtmp  = vt(1.00)
  λtmp  = -0.01
  cols  = fill(cm(0),length(Mode_Ind))
  for i in 1:length(Mode_Ind)
    cols[i] = cm(i-1)
  end

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
        jj = Int(i/histstep)
        Hist_Mode[jj,:,ik] = copy(z)
        Time[jj]           = t
    
        # Get field value at point x = hist_x,
        # corresponding to array index hist_i
        if (xhist)
          Histx[jj,ik] = Get_AsymptoticFieldx(hist_i,z,Vext,Y2,Y3)
        end  
      end

      if ifplot && mod(i,plotstep)==0

        # Remove previous plots
        for lo in ax3.get_lines()
          lo.remove()
        end
        for k in 1:length(Mode_Ind)
          ax3.plot(Time[1:jj],real.(Hist_Mode[1:jj,Mode_Ind[k],ik]),color=cm(k-1))
        end  

        # Remove previous plots
        for lo in ax4.get_lines()
          lo.remove()
        end
        ax4.plot(Time[1:jj],real.(Histx[1:jj,ik]),color=cm(ik-1))

        if (moveaxis)
          tmax = Time[jj]
          tmin = max(0.0,tmax-500.0)
          ax3.set_xlim([tmin,tmax])
          ax4.set_xlim([tmin,tmax])
        end  

        if (plotfield)
          # Remove previous plots
          for lo in ax5.get_lines()
            lo.remove()
          end
          fld12 = CenterManifold.GetAsymptoticField3(z,Vext,Y2,Y3)
          fld1  = fld12[1:ndof]
          ax5.plot(xg,real.(fld1),color=cm(0),linestyle="-", linewidth=1)
          ax5.plot(xg,imag.(fld1),color=cm(0),linestyle="--",linewidth=1)
          ax5.plot(xg,abs.(fld1), color=cm(1),linestyle="-", linewidth=3)
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


# Plot Peaks
if nsteps>0 && histplot
  h6          = figure(num=6,figsize=Grh.figsz1);
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
  ax6.plot(θA,Peak_Amp,linestyle="none",marker="o",markersize=Grh.mksz)
  ax6.set_xlabel(L"|θ|",fontsize=Grh.lafs)
  ax6.set_ylabel(L"A_{x}^{max}",fontsize=Grh.lafs)
end  


if (ifsave && nsteps>0)
  if ifresonant
    fname = "SL_pert_resonant_Parametric1.jld2"
  else
    fname = "SL_pert_nonresonant_Parametric1.jld2"
  end
  save(fname,"xg",xg,"Vext",Vext,"Y2",Y2,"Y3",Y3,"G1",G1,"G2",G2,"G3",G3,"δ",δ,"Time",Time,"θA",θA,"Peak_Amp",Peak_Amp,"Histx",Histx,"ω_nonlinear",ω_nonlinear,"Hist_Mode",Hist_Mode);
  println(fname*" saved.")
end 


println("Done.")












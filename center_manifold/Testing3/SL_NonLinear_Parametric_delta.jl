println("Non-Linear evolution for the Stuart Landau equations")

using Peaks
using Statistics
using Random
using JLD2

include("$JULIACOMMON/RK4.jl")
include("NLGinzburgLandau.jl")
include("OP_RK4.jl")

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
nsteps      = 3000000

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
δ4          = -Vector(1:10)*0.02
#δ4          = [-0.4]
nδ4         = length(δ4)
ncycles     = ones(Int64,nδ4)
# for i in 1:2
#   ncycles[i]  = 4
# end
# for i in 3:4
#   ncycles[i]  = 2
# end

Hist_Mode   = zeros(vt,nhist,m,nδ4)
Time        = zeros(Float64,nhist)
Peak_Amp    = zeros(Float64,nδ4)
ω_nonlinear = zeros(Float64,nδ4)
Mode_Ind    = [1]                         # Which mode to plot 


h3          = figure(num=3,figsize=Grh.figsz3);
ax3         = gca()
ax3.set_xlabel(L"time",fontsize=Grh.lafs)
ax3.set_ylabel(L"Z_{i}",fontsize=Grh.lafs)

TLast       = Tend - 500.0

if xhist
  hist_x,tmp= ForcingParams()             # Location of history point
  hist_i    = argmin(abs.(xg .- hist_x))  # Index of history point

  Histx     = zeros(vt,nhist,nδ4)

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
if (MaxOrd == 5)
  SL5(x)      = StuartLandau5(G1,G2,G3,G4,G5,x)
end

println("Press x to stop. Any other key to continue...")
xin = readline()
if xin == "x"
  nsteps = 0
end  


println("Starting Iterations")

for ik in 1:nδ4

  # Initialize z
  δ4_diff     = δ4[ik]

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
    j = i - (nsys+npert)
    if mod(j-1,2) == 0 
      z[i]    = δ4_diff
    else  
      z[i]    = z[i-1]'
    end
  end
  # Harmonic Forcing Amplitude
  for i in nsys+npert+p+1:m
    j = i - (nsys+npert+p)
    if mod(j-1,2) == 0 
      z[i]    = 0*(1.0 + 0.0im)
    else  
      z[i]    = z[i-1]'
    end
  end
  
  # Work Arrays
  zwork       = zeros(vt,m,5)

  cycles = ncycles[ik]
  for ic in 1:cycles

    # Start iterations
    t   = Inp.Dtype(0)    # Time
    
    for i in 1:nsteps
    
      t = t + dt;

      # Stuart Landau Evolution
      if (MaxOrd == 5)
        OP_RK4!(SL5,z,dt,zwork)
      else
        OP_RK4!(SL,z,dt,zwork)
      end

      # Set conjugation correctly
      # z[2] = z[1]'
    
      znorm = sqrt(abs(z'*z))
      # Print something  
      if verbose && mod(i,verbosestep)==0
        println("ik=$ik/$nδ4, ic=$ic/$cycles, Istep=$i, Time=$t, |z|=$znorm")
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
        ax3.plot(Time[1:jj],real.(Hist_Mode[1:jj,Mode_Ind,ik]),color=cm(ik-1))

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
          if MaxOrd == 5
            fld12 = CenterManifold.GetAsymptoticField(z,YM)
          else
            fld12 = CenterManifold.GetAsymptoticField3(z,Vext,Y2,Y3)
          end
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
    @printf("δ4: %.2f ; Amax: %.5f ; Ω: %.4e\n", δ4[ik], Peak_Amp[ik], ω_nonlinear[ik])
  end       # if nsteps>0 && histplot

end         # ik in 1:nδ4


# Plot Peaks
if nsteps>0 && histplot
  h6          = figure(num=6,figsize=Grh.figsz1);
  ax6         = gca()
  last_inds   = Time .> TLast
  Time2       = Time[last_inds]
  Hist2       = Histx[last_inds,:]
  Peak_Amp    = zeros(Float64,nδ4)
  for ik in 1:nδ4
    peak_ind        = argmaxima(real.(Hist2[:,ik]))
    peak_times      = Time2[peak_ind]
    maxamp          = real.(Hist2[peak_ind,ik])
    Peak_Amp[ik]    = mean(maxamp) 
    delta_times     = diff(peak_times)
    afreq           = 2.0*π./delta_times
    ω_nonlinear[ik] = mean(afreq)

    @printf("δ4: %.2f ; Amax: %.5f ; Ω: %.4e\n", δ4[ik],Peak_Amp[ik], ω_nonlinear[ik])
  end
  ax6.plot(abs.(δ4),Peak_Amp,linestyle="none",marker="o",markersize=Grh.mksz)
  ax6.set_xlabel(L"-δ4",fontsize=Grh.lafs)
  ax6.set_ylabel(L"A_{x}^{max}",fontsize=Grh.lafs)
end  


if (ifsave && nsteps>0)
  fname = "SL_diffusion_Parametric2.jld2"

  save(fname,"xg",xg,"Vext",Vext,"Y2",Y2,"Y3",Y3,"G1",G1,"G2",G2,"G3",G3,"δ",δ,"Time",Time,"δ4",δ4,"Peak_Amp",Peak_Amp,"Histx",Histx,"ω_nonlinear",ω_nonlinear,"Hist_Mode",Hist_Mode);
  println(fname*" saved.")
end 


println("Done.")












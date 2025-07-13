println("Time-Stepping Initialization")

using PolynomialBases
using PyPlot,PyCall
using LinearAlgebra
using IterativeSolvers

# Include the function files
include("sem_init_ref.jl")
include("custom_params.jl")
include("sem_main.jl")
include("$SRC/Dealias.jl")
include("$JULIACOMMON/GetEXT.jl")
include("$JULIACOMMON/GetBDF.jl")
include("$JULIACOMMON/GetBDF.jl")
include("RK4.jl")

fld         = zeros(VT,ndof,nflds)
dotfld      = zeros(VT,ndof,nflds)
fldlag      = zeros(VT,ndof,3,nflds)

Rhs         = zeros(VT,ndof,nflds)
Rhslag      = zeros(VT,ndof,2,nflds)

# Initial Condition
for j in 1:nflds

  # Gaussian
  for i in 1:ngauss
    amp       = Amp0[j]*ampgauss[i]
    fld[:,j]  = fld[:,j] .+ amp*vimultg.*QT*exp.(-0.5*((Geom.xm1[:] .- x0gauss[i])/σg).^2)
  end

  # Wavenumber
  for i in 1:nk0
    nk          = Normk0[i]
    amp         = ampk0[i,j]
    K0          = 2.0*π/(xe-xs)
    fld[:,j]    = fld[:,j] .+ amp*vimultg.*QT*sin.(K0*nk*Geom.xm1[:])
  end

  # Homogeneous and Noise
  fld[:,j]      = fld[:,j] .+ Off0[j] .+ σ0i[j]*(rand(ndof) .- 0.5)

end  


fldhist     = zeros(VT,npts,nsurf_save,nflds)
Thist       = zeros(VT,nsurf_save)

bdf         = zeros(Float64,4)
bdfto2      = zeros(Float64,4)
ext         = zeros(Float64,3)

cm          = get_cmap("tab10");
rgba0       = cm(0); 
rgba1       = cm(1); 
rgba2       = cm(2);

time        = range(0.,step=dt,length=nsteps);

h2          = figure(num=2)
ax2         = h2.subplots()
ax2.set_xlim(xs,xe)
ax2.set_ylim(-1.5,6.0)
ax2.set_xlabel(L"x", fontsize=lafs)
ax2.set_ylabel(L"u,v", fontsize=lafs)
MoveFigure(h2,1250,10)

t           = 0.

framecount  = 0


Thist[1] = t
for i in 1:nflds
  fldhist[:,1,i] = Q*fld[:,i]
end  

# Plot Initial Conditions
pl = Array{Any}(undef,nflds)
if initplot
  for j in 1:nflds
    if (plotfldi[j])
      if j == 1
        col = cmapV
      elseif j == 2
        col = cmapU
      else
        col = cmap(j-1)
      end  
      pl[j] = ax2.plot(Geom.xm1[:],Q*fld[:,j],color=col);
#      pl2   = ax2.plot(Geom.xm1[:],Q*fld[:,2],color=rgba1);
    end
  end  
  if ifphplot
    scat = ax1.plot(fld[:,1],fld[:,2],color="black") 
  end

  if (ifsaveframe)
    framecount = framecount + 1
    if (ifphplot)
      fname   = @sprintf "./plots/phase/phase_%06i" framecount
      h1.savefig(fname)
    end

    if (iffldplot)
      fname   = @sprintf "./plots/fields/fields_%06i" framecount
      h2.savefig(fname)
    end
  end  

  println("Press x to stop. Any other key to continue")
  StepperStart = readline()
  if StepperStart !="x"
    # Remove Plots  
    for j in 1:nflds
      if (plotfldi[j])
        pl[j][1].remove();
      end
    end  
    if (ifphplot)
      scat[1].remove()
    end 
    # Dynamically moving Null-clines
    if (ifdynnull)
      for i in 1:length(PlotContainers)
        if isassigned(PlotContainers,i)
          p = PlotContainers[i]
          p[1].remove()
          PlotContainers[i] = []
        end  
      end
    end

    # Start Stepper
    include("time_stepper_multiple_FitzHughNagumo.jl")
  end       # StepperStart != "x"

end

# # Plot Dynamically moving Null-clines
# if (ifdynnull)
#   for i in 1:length(PlotContainers)
#     if isassigned(PlotContainers,i)
#       p = PlotContainers[i]
#       p[1].remove()
#       PlotContainers[i] = []
#     end  
#   end
# end

pause(0.01)






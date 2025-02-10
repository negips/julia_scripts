println("Time-Stepping interface")

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

include("time_stepper_multiple_init.jl")

X = Geom.xm1[:];
X[end] = X[1]
QTX = QT*(X.*vimult)

# No dynamic phase here
ifdynplot         = true
ifplot            = iffldplot || ifphplot

Vol   = sum(Bg)
A_sen = Asen*Vol
γ     = (2.0)/λnorm

γhist       = zeros(VT,nsurf_save)
Abarhist    = zeros(VT,nsurf_save)
γhist[1]    = γ
Abarhist[1] = (Bg'*fld[:,2])/A_sen - Aeq


for i in 1:nsteps
  global fld,fldlag,Rhs,Rhslag,dotfld
  global t
  global pl,pl2,scat,λpl
  global PlotContainers
  global framecount
  global γ

  t = t + dt;

  A_tot     = Bg'*fld[:,2]
  abar      = A_tot/A_sen - Aeq

  γ         = RK4!(λdot1,abar,γ,dt)

  if verbosestep>0 && mod(i,verbosestep)==0
    println("Step: $i/$nsteps, Time: $t, Abar = $(abar); γ: $(γ)")
  end

  GetBDF!(bdf,3)
  GetEXT!(ext,3)
    
  if i==1
    GetBDF!(bdf,1)
    GetEXT!(ext,1)
  elseif i==2
    GetBDF!(bdf,2)
    GetEXT!(ext,2)
  end

  θpar      = (θ0 + (γ/1.0)*dθ)*π/180.0
  #Ω         = 0.05
  #θpar      = (θ0 - (A_tot/A_eq)*dθ)*pi/180.0
  # θpar      = 0.0 # (θ0 + dθ*sin(2*π*Ω*t))*pi/180.0         # for G
  λpar      = 0.0 # (λ0 - dλ*sin(2*π*Ω*t))                  # for F
  dotfld    = Flow(fld[:,1],fld[:,2],θpar,λpar)

  for j in 1:nflds
    rhs           = dotfld[:,j] .- Filg*fld[:,j];
    rhs1          = ext[1]*rhs + ext[2]*Rhslag[:,1,j] + ext[3]*Rhslag[:,2,j];

    Rhslag[:,2,j] = copy(Rhslag[:,1,j]);
    Rhslag[:,1,j] = copy(rhs);

    bdlag         = 1. /dt*(bdf[2]*fld[:,j] + bdf[3]*fldlag[:,1,j] + bdf[4]*fldlag[:,2,j]); 
    Rhs[:,j]      = Bg.*(rhs1 .- bdlag);

    fldlag[:,2,j] = copy(fldlag[:,1,j]);
    fldlag[:,1,j] = copy(fld[:,j]);

#   Noise    
    Σ             = σall[j]*(rand(ndof) .- 0.5)
    Rhs[:,j]      = Rhs[:,j] .+ Bg.*Σ
  end

  for j in 1:nflds
    M         = bdf[1]/dt*diagm(Bg) .- γall[j]*Lg;
    a         = gmres(M,Rhs[:,j],abstol=1.0e-10,verbose=false)
    fld[:,j]  = copy(a)
  end

  # Save for surface plot  
  if surf_save>0 && mod(i,surf_save)==0
    k = Int(i/surf_save) + 1
    Thist[k]  = t
    for j in 1:nflds
      fldhist[:,k,j] = Q*fld[:,j]
    end

    γhist[k]   = γ
    Abarhist[k]= abar
  end  

  if ifplot && mod(i,plotupd)==0
    # Remove old plots      
    if (i>plotupd)
       if (iffldplot)
         for j in 1:nflds
          if (plotfldi[j])
            pl[j][1].remove()
          end
        end  
       end  
       if (ifphplot)
         scat[1].remove()
       end
       if ifdynplot
         λpl[1].remove()
       end

       if (ifdynnull)
#         PlotContainers[1][1].remove()
#         PlotContainers[2][1].remove()
#
#         PlotContainers[3][1].remove()
#         PlotContainers[4][1].remove()

         PlotContainers[5][1].remove()
         PlotContainers[6][1].remove()
#         PlotContainers[7][1].remove()
       end  
      
    end

    # Add updated plots      
    if (iffldplot)
      for j in 1:nflds
        if (plotfldi[j])
          pl[j] = ax2.plot(Geom.xm1[:],Q*fld[:,j],color=cm(j-1));
        end
      end  
    end
    # Phase plot    
    if ifphplot
      scat = ax1.plot(Intpg*fld[:,1],Intpg*fld[:,2],color="black",linewidth=2) 
    end
   
    # Dynamic plot
    if ifdynplot
      # λpl =ax2.plot(Geom.xm1[:],Q*λpar,color=cm(4-1));
      λpl = ax4.plot(abar,γ,color=cm(4-1),linestyle="none",marker="o",markersize=8);
    end

    # Dynamic null-clines
    if ifdynnull
      PlotContainers[5] = ax1.plot(ft(λpar),yin,linestyle="--",linewidth=2,color=cm(0));
      PlotContainers[6] = ax1.plot(gt(θpar),yin,linestyle="--",linewidth=2,color=cm(1));
    end  

    # Saving frames    
    if (ifsaveframe)
      framecount = framecount + 1
      if (ifphplot)
        fname2   = @sprintf "./plots/phase/phase_%06i" framecount
        h1.savefig(fname2)
      end

      if (iffldplot)
        fname2   = @sprintf "./plots/fields/fields_%06i" framecount
        h2.savefig(fname2)
      end
    end  

    pause(0.001)
  end  

end


t2d   = ones(npts)*Thist'
x2d   = (Geom.xm1[:])*ones(nsurf_save)'

cm2   = get_cmap("binary");
h3    = figure(num=3,figsize=[8.0,8.0])
pcm   = pcolormesh(x2d,t2d,fldhist[:,:,2])
pcm.set_cmap(cm2)
ax3   = h3.gca()
ax3.invert_yaxis()
#cb    = colorbar(orientation="vertical")
if (ifsavext)
  fname3 = @sprintf "./plots/spacetime"
  h3.savefig(fname3)
  println("Saved Figure "*fname3)
end  
#surf(t2d,x2d,fldhist[:,:,2],cmap=cm2,edgecolor="none")
#ax3.elev = 94.0
#ax3.azim = 0.0
#ax3.roll = 0.0
#draw()

h5,(ax5,ax6) = subplots(1,2,sharey=true,figsize=[8.0,8.0])
ax5.set_position([0.125, 0.10, 0.55, 0.8])
ax6.set_position([0.700, 0.10, 0.20, 0.8])
sca(ax5)
pcm   = pcolormesh(x2d,t2d,fldhist[:,:,2])
pcm.set_cmap(cm2)
ax5.set_ylabel("t",fontsize=lafs)
ax5.set_xlabel("x",fontsize=lafs)
ax5.invert_yaxis()

ax6.plot(γhist,Thist,linewidth=2,color=cm(3))
ax6.set_xlabel("λ",fontsize=lafs)
if (ifsavext)
  fname4   = @sprintf "./plots/spacetime2"
  h5.savefig(fname4)
  println("Saved Figure "*fname4)
end  








println("Time-Stepping interface")

using PolynomialBases
using PyPlot,PyCall
using LinearAlgebra
using IterativeSolvers

#close("all")

# Include the function files
include("sem_main.jl")
#include("Meinhardt.jl")
include("Dealias.jl")
include("../GetEXT.jl")
include("../GetBDF.jl")

include("time_stepper_multiple_init.jl")

X = Geom.xm1[:];
X[end] = X[1]
QTX = QT*(X.*vimult)

# No dynamic phase here
ifdynplot         = true
ifplot            = iffldplot || ifphplot

for i in 1:nsteps
  global fld,fldlag,Rhs,Rhslag,dotfld
  global t
  global pl,pl2,scat,λpl
  global PlotContainers
  global framecount

  t = t + dt;

  if verbosestep>0 && mod(i,verbosestep)==0
    println("Step: $i, Time: $t")
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

  dotfld    = Flow(fld[:,1],fld[:,2],fld[:,3])
#  λpar      = 0.25*copy(fld[:,3])
  λpar      = (fld[:,2]./eq_a).^2 .- 1.5*(fld[:,3]/eq_λ)
#  λpar      = ((fld[:,2]/eq_a) .- 1.1*(fld[:,3]./eq_λ))
  dotfld2   = Flow(fld[:,1],fld[:,2],1.0*λpar)
  dotfld[:,1:2] = dotfld2[:,1:2]

  for j in 1:nflds
    rhs           =  dotfld[:,j] .- Filg*fld[:,j];
    rhs1          =  ext[1]*rhs + ext[2]*Rhslag[:,1,j] + ext[3]*Rhslag[:,2,j];

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
  end  

  if ifplot && mod(i,plotupd)==0
#   Remove old plots      
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
         PlotContainers[7][1].remove()
       end  
      
    end

#   Add updated plots      
    if (iffldplot)
      for j in 1:nflds
        if (plotfldi[j])
          pl[j] = ax2.plot(Geom.xm1[:],Q*fld[:,j],color=cm(j-1));
        end
      end  
    end
#   Phase plot    
    if ifphplot
      scat = ax1.plot(fld[:,1],fld[:,2],color="black") 
    end
    
    if ifdynplot
      λpl =ax2.plot(Geom.xm1[:],Q*λpar,color=cm(4-1));
    end

    if ifdynnull
#      λmin = minimum(λpar)
#      PlotContainers[1] = ax1.plot(ft(λmin),yin,linestyle=":",linewidth=2,color=cm(0));
#      PlotContainers[2] = ax1.plot(gt(λmin),yin,linestyle=":",linewidth=2,color=cm(1));

#      λmax = maximum(λpar)
#      PlotContainers[3] = ax1.plot(ft(λmax),yin,linestyle="--",linewidth=2,color=cm(0));
#      PlotContainers[4] = ax1.plot(gt(λmax),yin,linestyle="--",linewidth=2,color=cm(1));

      index = argmax(fld[:,3])
#      index = argmin(abs.(QTX .- x0gauss[1]))      
      λmax  = λpar[index]
      PlotContainers[5] = ax1.plot(ft(λmax),yin,linestyle="--",linewidth=2,color=cm(0));
      PlotContainers[6] = ax1.plot(gt(λmax),yin,linestyle="--",linewidth=2,color=cm(1));

      bmax  = fld[index,1]
      amax  = fld[index,2]
      PlotContainers[7] = ax1.plot(bmax,amax,linestyle="none",marker="s",markersize=6,color=cm(2));
    
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
h3    = figure(num=3)
pcm   = pcolormesh(x2d,t2d,fldhist[:,:,2])
pcm.set_cmap(cm2)
ax3   = h3.gca()
ax3.invert_yaxis()
cb    = colorbar(orientation="vertical")

#surf(t2d,x2d,fldhist[:,:,2],cmap=cm2,edgecolor="none")
#ax3.elev = 94.0
#ax3.azim = 0.0
#ax3.roll = 0.0
#draw()









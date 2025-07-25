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
ifdynplot         = false
ifplot            = iffldplot || ifphplot

Vol  = sum(Bg)

for i in 1:nsteps
  global fld,fldlag,Rhs,Rhslag,dotfld
  global t
  global pl,pl2,scat,λpl
  global PlotContainers
  global framecount

  t = t + dt;

  A_tot     = Bg'*fld[:,2]/Vol

  if verbosestep>0 && mod(i,verbosestep)==0
    println("Step: $i/$nsteps, Time: $t")
  end

  GetBDF!(bdf,3)
  GetEXT!(ext,3)
    
  if i==1
    GetBDF!(bdf,1)
    GetBDFTO2!(bdfto2,1)
    GetEXT!(ext,1)
  elseif i==2
    GetBDF!(bdf,2)
    GetBDFTO2!(bdfto2,2)
    GetEXT!(ext,2)
  end

  dotfld    = Flow(fld[:,1],fld[:,2])

  for j in 1:nflds
    rhs           = dotfld[:,j] .- Filg*fld[:,j] .+ Cg*fld[:,j]
    rhs1          = ext[1]*rhs + ext[2]*Rhslag[:,1,j] + ext[3]*Rhslag[:,2,j]

    Rhslag[:,2,j] = copy(Rhslag[:,1,j])
    Rhslag[:,1,j] = copy(rhs)

    bdlag         = 1.0/dt*(bdf[2]*fld[:,j] + bdf[3]*fldlag[:,1,j] + bdf[4]*fldlag[:,2,j]) 
    Rhs[:,j]      = Bg.*(rhs1 .- bdlag)

    bdlag2        = 1.0/(dt^2)*(bdfto2[2]*fld[:,j] + bdfto2[3]*fldlag[:,1,j] + bdfto2[4]*fldlag[:,2,j])
    Rhs[:,j]      = Bg.*(rhs1 .- bdlag2)
 
    fldlag[:,3,j] = copy(fldlag[:,2,j])
    fldlag[:,2,j] = copy(fldlag[:,1,j])
    fldlag[:,1,j] = copy(fld[:,j])

    # Noise    
    Σ             = σall[j]*(rand(ndof) .- 0.5)
    Rhs[:,j]      = Rhs[:,j] .+ Bg.*Σ
  end

  for j in 1:nflds
    M         = (bdf[1]/dt + bdfto2[1]/dt/dt)*diagm(Bg) .- γall[j]*Lg;
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
    # Remove old plots      
    if (i>plotupd)
       if (iffldplot)
         for j in 1:nflds
          if (plotfldi[j])
            pl[j][1].remove()
          end
        end
        pl[nflds+1][1].remove()
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
      pl[nflds+1] = ax2.plot(Geom.xm1[:],Q*(fld[:,1] .+ fld[:,2]),color=cm(nflds));
    end
    # Phase plot    
    if ifphplot
      scat = ax1.plot(Intpg*fld[:,1],Intpg*fld[:,2],color="gray",linewidth=2) 
    end
   
    # Dynamic plot
    if ifdynplot
      λpl =ax2.plot(Geom.xm1[:],Q*λpar,color=cm(4-1));
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

if (ftype == 1)
  vmin = -1.0
  vmax =  1.0
elseif (ftype == 2)
  vmin =  0.0
  vmax =  1.0
elseif (ftype == 3)
  vmin = -1.0
  vmax =  1.0
end  


#cm2   = get_cmap("binary");
cm2   = get_cmap("seismic") # coolwarm, bwr, seismic
h3    = figure(num=3,figsize=[5.0,8.0])
pcm   = pcolormesh(x2d,t2d,fldhist[:,:,2],vmin=vmin,vmax=vmax)
pcm.set_cmap(cm2)
ax3   = h3.gca()
ax3.invert_yaxis()
# cb    = colorbar(orientation="vertical")
if (ifsavext)
  fname3 = @sprintf "./plots/spacetime"
  h3.savefig(fname3)
end  

#surf(t2d,x2d,fldhist[:,:,2],cmap=cm2,edgecolor="none")
#ax3.elev = 94.0
#ax3.azim = 0.0
#ax3.roll = 0.0
#draw()









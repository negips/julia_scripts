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


include("time_stepping_multiple_init.jl")


# Dynamic variables
#------------------------------ 
vol         = sum(Bg)
λ           = 0.1
θ           = sum(Bg.*fld[:,2])/vol
θ0          = 1.0               # Equilibrium mean Activation
τ           = 10.0              # Time scale of change
τi          = 1.0/τ             # Inverse τ

λlag        = zeros(VT,ndof,3)
λRhslag     = zeros(VT,ndof,3)
#------------------------------ 

#h4  = figure(num=4)
#ax4 = h4.subplots()
#λpl0 = ax4.plot(θ,λ,markersize=6,marker="o",color="black",fillstyle="none")
#ax4.set_xlabel(L"θ", fontsize=lafs)
#ax4.set_ylabel(L"λ", fontsize=lafs)
#ax4.set_xlim(-0.2,2.0)
#ax4.set_ylim(0.0,5.0)

#MoveFigure(h4,600,500)

if !ifdynplot
  close(h4)
end

for i in 1:nsteps
  global fld,fldlag,Rhs,Rhslag,dotfld
  global t
  global pl,pl2,scat
  global λ,λlag,λRhslag,θ
  global h4,ax4,λpl

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

# Solve for the dynamic variable λ
#----------------------------------------
  active          = fld[:,2] .> 0.0
  θ               = 3.0*sum(Bg.*active.*fld[:,2])/vol
  rhsλ            = λdot(θ,λ) 
  rhs1λ           = ext[1]*rhsλ + ext[2]*λRhslag[1] + ext[3]*λRhslag[2];

  λRhslag[2]      = copy(λRhslag[1]);
  λRhslag[1]      = copy(rhsλ);

  bdλ             = 1. /dt*(bdf[2]*λ + bdf[3]*λlag[1] + bdf[4]*λlag[2]); 
  λlag[2]         = copy(λlag[1]);
  λlag[1]         = copy(λ);
 
  Rλ              = rhs1λ - bdλ;
  λ               = Rλ*dt/bdf[1]
#----------------------------------------   

  dotfld = 1.0(λ*Flow1(fld[:,1],fld[:,2]) .+ (1.0 - λ)*Flow2(fld[:,1],fld[:,2]))

# Build Rhs  
  for j in 1:nflds
    rhs           =  dotfld[:,j] .- Filg*fld[:,j];
    rhs1          =  ext[1]*rhs .+ ext[2]*Rhslag[:,1,j] .+ ext[3]*Rhslag[:,2,j];

    Rhslag[:,2,j] = copy(Rhslag[:,1,j]);
    Rhslag[:,1,j] = copy(rhs);

    bdlag         = 1. /dt*(bdf[2]*fld[:,j] .+ bdf[3]*fldlag[:,1,j] .+ bdf[4]*fldlag[:,2,j]); 
    Rhs[:,j]      = Bg.*(rhs1 .- bdlag);

    fldlag[:,2,j] = copy(fldlag[:,1,j]);
    fldlag[:,1,j] = copy(fld[:,j]);
    
#   Noise    
    Σ             = σall[j]*(rand(ndof) .- 0.5)
    Rhs[:,j]      = Rhs[:,j] .+ Bg.*Σ
  end  

# Solve  
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

# Dynamic plot update
  if ifplot && mod(i,plotupd)==0
    if (i>plotupd)
#      Field plots         
       if (iffldplot)
         pl[1].remove()
         pl2[1].remove()
       end
#      Phase Plots         
       if (ifphplot)
         scat[1].remove()
       end
#      Dynamic phase (λ)         
       if ifdynplot
         λpl[1].remove()
       end
    end   
#   Updated plots      
    if (iffldplot)
      pl  = ax2.plot(Geom.xm1[:],Q*fld[:,1],color=rgba0)
      pl2 = ax2.plot(Geom.xm1[:],Q*fld[:,2],color=rgba1)
    end  
    if ifphplot
      scat = ax1.plot(fld[:,1],fld[:,2],color="black") 
    end
    if ifdynplot
      λpl = ax4.plot(θ,λ,markersize=6,marker="o",color="black")
    end
    draw()
    pause(0.001)
  end  

end


cm2   = get_cmap("binary");
t2d   = ones(npts)*Thist'
x2d   = (Geom.xm1[:])*ones(nsurf_save)'

h3  = figure(num=3)
pcm = pcolormesh(x2d,t2d,fldhist[:,:,2])
pcm.set_cmap(cm2)
ax3 = h3.gca()
ax3.invert_yaxis()
cb  = colorbar()

#surf(t2d,x2d,fldhist[:,:,2],cmap=cm2,edgecolor="none")
#ax3.elev = 90.0
#ax3.azim = 0.0
#draw()













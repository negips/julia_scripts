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

agauss      = exp.(-((Geom.xm1[:] .- x0)/σg).^2)# .*(sign.(Geom.xm1[:] .- x0))
k0          = 0.5
asin        = sin.(k0*Geom.xm1[:])

ainit       = vimultg.*(QT*agauss)

nflds = 2                                 # No of fields
fld     = zeros(VT,ndof,nflds)
dotfld  = zeros(VT,ndof,nflds)
fldlag  = zeros(VT,ndof,2,nflds)

Rhs     = zeros(VT,ndof,nflds)
Rhslag  = zeros(VT,ndof,2,nflds)

fld[:,1] = ampB0*ainit .+ B0Off
fld[:,2] = ampA0*ainit .+ A0Off

fldhist  = zeros(VT,nsurf_save,npts,nflds)
Thist    = zeros(VT,nsurf_save)

bdf = zeros(Float64,4)
ext = zeros(Float64,3)

cm    = get_cmap("tab10");
rgba0 = cm(0); 
rgba1 = cm(1); 
rgba2 = cm(2); 

time = range(0.,step=dt,length=nsteps);

h2  = figure(num=2)
ax2 = h2.subplots()

t = 0.

Thist[1] = t
for i in 1:nflds
  fldhist[1,:,i] = Q*fld[:,i]
end  

for i in 1:nsteps
  global fld,fldlag,Rhs,Rhslag,dotfld
  global t
  global pl,pl2,scat

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

  if i == nstep_switch
    println("Switching Functions.")
  end  
  if i<nstep_switch
    dotfld = Flow1(fld[:,1],fld[:,2])
  else
    dotfld = Flow2(fld[:,1],fld[:,2])
  end

  for j in 1:nflds
    rhs           =  dotfld[:,j] .- Filg*fld[:,j];
    rhs1          =  ext[1]*rhs + ext[2]*Rhslag[:,1,j] + ext[3]*Rhslag[:,2,j];

    Rhslag[:,2,j] = copy(Rhslag[:,1,j]);
    Rhslag[:,1,j] = copy(rhs);

    bdlag         = 1. /dt*(bdf[2]*fld[:,j] + bdf[3]*fldlag[:,1,j] + bdf[4]*fldlag[:,2,j]); 
    Rhs[:,j]      = Bg.*(rhs1 .- bdlag);

    fldlag[:,2,j] = copy(fldlag[:,1,j]);
    fldlag[:,1,j] = copy(fld[:,j]);
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
      fldhist[k,:,j] = Q*fld[:,j]
    end
  end

  if plotupd > 0
    if mod(i,plotupd)==0
      if (i>plotupd)
         pl[1].remove()
         pl2[1].remove()
         if (ifphplot)
           scat[1].remove()
         end  
      end   
      pl = ax2.plot(Geom.xm1[:],Q*fld[:,1],color=rgba0);
#      pl = plot(Geom.xm1[:],Q*(fld[:,1]),color=rgba1);
     
      pl2 = ax2.plot(Geom.xm1[:],Q*fld[:,2],color=rgba1);

      k     = argmax(abs.(fld[:,2]))
      xscat = fld[k,1]
      yscat = fld[k,2]

      if ifphplot
#        scat  = ax1.plot(xscat,yscat,marker="o",color="black")
        scat = ax1.plot(fld[:,1],fld[:,2],color="black") 
      end
      draw()
      pause(0.01)
    end
  end  

end


cm2   = get_cmap("Spectral");
t2d   = Thist*ones(1,npts)
x2d   = ones(nsurf_save)*(Geom.xm1[:]')

h3=figure(num=3)
surf(t2d,x2d,fldhist[:,:,2],cmap=cm2,edgecolor="none")
ax3 = h3.gca()
ax3.elev = 90.0
ax3.azim = 0.0
draw()













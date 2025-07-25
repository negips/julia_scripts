println("Main interface for 1D SEM")

using PolynomialBases
using PyPlot,PyCall
using LinearAlgebra
using IterativeSolvers
using Roots

close("all")

const SRC = "/home/prabal/workstation/git/julia/OIST/Custom"

# Include the function files
include("sem_main.jl")
include("FitzhughNagumo.jl")
include("Dealias.jl")
include("$JULIACOMMON/GetEXT.jl")
include("$JULIACOMMON/GetBDF.jl")
include("$JULIACOMMON/MoveFigure.jl")

FHN(u,v) = FitzhughNagumo(u,v,a,b)
include("plot_fitzhughnagumo_nullclines.jl")

aa    = exp.(-((Geom.xm1[:] .- x0)/σ).^2)
ainit = vimultg.*(QT*aa)

nflds = 2                                 # No of fields
fld     = zeros(VT,ndof,nflds)
dotfld  = zeros(VT,ndof,nflds)
fldlag  = zeros(VT,ndof,2,nflds)

Rhs     = zeros(VT,ndof,nflds)
Rhslag  = zeros(VT,ndof,2,nflds)

#Qfld    = zeros(VT,npts)
#Qcfld   = zeros(VT,npts)
#Qconv   = zeros(VT,npts)
#convect = zeros(VT,ndof)

fld[:,1] = ampU0*copy(ainit) .+ U0Off
fld[:,2] = ampV0*copy(ainit) .+ V0Off

#blag  = zeros(Float64,ndof,3);
#Rblag = zeros(Float64,ndof,2);

bdf   = zeros(Float64,4)
ext   = zeros(Float64,3)

cm    = get_cmap("tab10");
rgba0 = cm(0); 
rgba1 = cm(1); 
rgba2 = cm(2); 

time  = range(0.,step=dt,length=nsteps);

h2    = figure(num=2)
ax2   = h2.subplots()

println("Press any key to continue")
xinn  = readline()

t   = 0.
for i in 1:nsteps
  global fld,fldlag,Rhs,Rhslag,dotfld
  global t
  global pl,pl2
  global ax1,ax2
  global scat

  t = t + dt;

  GetBDF!(bdf,3)
  GetEXT!(ext,3)
    
  if i==1
    GetBDF!(bdf,1)
    GetEXT!(ext,1)
  elseif i==2
    GetBDF!(bdf,2)
    GetEXT!(ext,2)
  end

  dotfld = FHN(fld[:,1],fld[:,2])
  for j in 1:nflds
    # Qfld          = Q*fld[:,j]                        # Convected field  
    # Qcfld         = cfall[j]*Q*fld[:,1] .- ζall[j]    # Convecting field (is always u)
    # Dealias!(Qconv,Qcfld,Qfld,Geom)
    # convect       = (QT*Qconv)./Bg

    if j == 1
      dotfld[:,j] = (1.0/ϵ)*dotfld[:,j]
    end  

    rhs           =  dotfld[:,j] .- Filg*fld[:,j] # .- convect;
    rhs1          =  ext[1]*rhs + ext[2]*Rhslag[:,1,j] + ext[3]*Rhslag[:,2,j];

    Rhslag[:,2,j] = copy(Rhslag[:,1,j]);
    Rhslag[:,1,j] = copy(rhs);

    bdlag         = 1. /dt*(bdf[2]*fld[:,j] + bdf[3]*fldlag[:,1,j] + bdf[4]*fldlag[:,2,j]); 
    Rhs[:,j]      = Bg.*(rhs1 .- bdlag);

    fldlag[:,2,j] = copy(fldlag[:,1,j]);
    fldlag[:,1,j] = copy(fld[:,j]);
  end  

  for j in 1:nflds
    M        = bdf[1]/dt*diagm(Bg) .- γall[j]*Lg;
    local a  = gmres(M,Rhs[:,j],abstol=1.0e-10,verbose=false)
    fld[:,j] = copy(a)
  end 

  if plotupd > 0
    if mod(i,plotupd)==0
      if (i>plotupd)
        pl[1].remove();
        pl2[1].remove();
        if (ifphplot)
          scat[1].remove()
        end
      end

      pl = ax2.plot(Geom.xm1[:],Q*fld[:,1],color=rgba0);
#      pl = plot(Geom.xm1[:],Q*(fld[:,1]),color=rgba1);
      pl2 = ax2.plot(Geom.xm1[:],Q*fld[:,2],color=rgba1);

      if ifphplot
        scat = ax1.plot(Intpg*fld[:,1],Intpg*fld[:,2],color="black") 
      end
      pause(0.001)
    end     # mod(i,plotupd)
  end       # plotupd > 0 
end         # i in 1:nsteps















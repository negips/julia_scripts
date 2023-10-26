println("Main interface for 1D SEM")

using PolynomialBases
using PyPlot,PyCall
using LinearAlgebra
using IterativeSolvers

close("all")

# Include the function files
include("sem_main.jl")
include("Meinhardt.jl")
include("../GetEXT.jl")
include("../GetBDF.jl")


Mein(a,b) = Meinhardt_1987_2_branching(a,b,R)

x0 = 10.0
σ  = 1.0
aa = 1.0*exp.(-((Geom.xm1[:] .- x0)/σ).^2)
ainit = QT*(aa./vmult)
a     = copy(ainit);
alag  = zeros(Float64,ndof,3);
Ralag = zeros(Float64,ndof,2);

b     = 0.1*a
blag  = zeros(Float64,ndof,3);
Rblag = zeros(Float64,ndof,2);

#xall  = zeros(Float64,r);
#for i in 1:nel
#   global xall   
#   j1 = (i-1)*N + 1;
#   j2 = j1+lx1-1;
#   xall[j1:j2] = Geom.xm1[:,i];
#end   

bdf = zeros(Float64,4)
ext = zeros(Float64,3)

cm    = get_cmap("tab10");
rgba0 = cm(0); 
rgba1 = cm(1); 
rgba2 = cm(2); 

dt = 0.001;
plotupd = 100;

nsteps = 100000;

time = range(0.,step=dt,length=nsteps);

t = 0.
for i in 1:nsteps
  global a,b,alag,Ralag,blag,Rblag
  global t
  global pl,pl2

  t = t + dt;

  GetBDF!(bdf,3)
  GetEXT!(ext,3)
    
  if i==1
    GetBDF!(bdf,1)
    GetEXT!(ext,1)
#    bdf = bdf1;
#    ext = ex0;
  elseif i==2
    GetBDF!(bdf,2)
    GetEXT!(ext,2)
   
#    bdf = bdf2;
#    ext = ex1;
#  elseif i==3 
#    bdf = bdf3;
#    ext = ex2;
  end

# Extrapolate rhs
  adot,bdot     = Mein(a,b) 
  rhs     = adot .- Filg*a;
  rhsa    = ext[1]*rhs + ext[2]*Ralag[:,1] + ext[3]*Ralag[:,2];

  Ralag[:,2] = Ralag[:,1];
  Ralag[:,1] = rhs;

  bdlag    = 1. /dt*(bdf[2]*a + bdf[3]*alag[:,1] + bdf[4]*alag[:,2]); 
  Rhsa     = Bg.*(rhsa - bdlag);

  rhs     = bdot .- Filg*b;
  rhsb    = ext[1]*rhs + ext[2]*Rblag[:,1] + ext[3]*Rblag[:,2];

  Rblag[:,2] = Rblag[:,1];
  Rblag[:,1] = rhs;


  bdlag    = 1. /dt*(bdf[2]*b + bdf[3]*blag[:,1] + bdf[4]*blag[:,2]); 
  Rhsb     = Bg.*(rhsb - bdlag);

  alag[:,3] = alag[:,2];
  alag[:,2] = alag[:,1];
  alag[:,1] = a;

  blag[:,3] = blag[:,2];
  blag[:,2] = blag[:,1];
  blag[:,1] = b;

  Ma        = bdf[1]/dt*diagm(Bg) .- γa*Lg;
#  a         = Ma\Rhsa
   gmres!(a,Ma,Rhsa,abstol=1.0e-10,verbose=false)

  Mb        = bdf[1]/dt*diagm(Bg) .- γb*Lg;
  gmres!(b,Mb,Rhsb,abstol=1.0e-10)
#  b         = Mb\Rhsb
 
  if mod(i,plotupd)==0
    if (i>plotupd)
       pl[1].remove();
       pl2[1].remove();
    end   
  
    pl = plot(Geom.xm1[:],Q*a,color=rgba0);
    pl2 = plot(Geom.xm1[:],Q*b,color=rgba1);

    pause(0.001)
  end  

end
















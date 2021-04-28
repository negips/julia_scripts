println("Main interface for 1D SEM")

using PolynomialBases
using PyPlot,PyCall
using LinearAlgebra


close("all")

# Include the function files
include("sem_main.jl")


# Temporal discretization
bdf1  = [ 1.  -1.  0.  0.]/1.;
bdf2  = [ 3.  -4.  1.  0.]/2.;
bdf3  = [11. -18.  9. -2.]/6.;
ex0   = [0. 1.  0. 0.];
ex1   = [0. 2. -1. 0.];
ex2   = [0. 3. -3. 1.];

r,c   = size(L);

Vinit = rand(Float64,r);
V     = copy(Vinit);
Vlag  = zeros(Float64,r,3);
Rlag  = zeros(Float64,r,2);

xall  = zeros(Float64,r);
for i in 1:nel
   global xall   
   j1 = (i-1)*N + 1;
   j2 = j1+lx1-1;
   xall[j1:j2] = Geom.xm1[:,i];
end   


cm    = get_cmap("tab10");
rgba0 = cm(0); 
rgba1 = cm(1); 
rgba2 = cm(2); 

dt = 0.01;
plotupd = 10;

nsteps = 10000;

time = range(0.,step=dt,length=nsteps);

t = 0.
for i in 1:nsteps
  global V,Vlag,Rlag
  global t
  global pl

  t = t + dt;

  bdf = bdf3;
  ext = ex2;
     
  if i==1
    bdf = bdf1;
    ext = ex0;
  elseif i==2  
    bdf = bdf2;
    ext = ex1;
#  elseif i==3 
#    bdf = bdf3;
#    ext = ex2;
  end

  A = L;
  A[1,:] = zeros(Float64,1,r);
  A[1,1] = 1.

# Extrapolate A*x
  rhs     = A*V;
  rhs1    = ext[2]*rhs + ext[3]*Rlag[:,1] + ext[4]*Rlag[:,2];

  Rlag[:,2] = Rlag[:,1];
  Rlag[:,1] = rhs;

  bdlag    = 1. /dt*(bdf[2]*V + bdf[3]*Vlag[:,1] + bdf[4]*Vlag[:,2]); 
  Rhs      = rhs1 - bdlag;

  Vlag[:,3] = Vlag[:,2];
  Vlag[:,2] = Vlag[:,1];
  Vlag[:,1] = V;

  Rhs       = Rhs*dt/bdf[1];
  Rhs[1]    = 1.0*sin(t);

  V         = Rhs;

  if mod(i,plotupd)==0
    if (i>plotupd)
       pl[1].remove();
    end   
  
    pl = plot(xall,V,color=rgba0);
 
    pause(0.001)
  end  

end
















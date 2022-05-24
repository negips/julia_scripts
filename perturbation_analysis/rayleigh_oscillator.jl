println("Rayleigh Oscillator")

import Pkg
#Pkg.add("PyPlot")
#Pkg.add("Blink")

using Blink
using PyPlot,Colors,PyCall
using LinearAlgebra
using Random

lafs = 16;

# Temporal discretization
bdf1  = [ 1.  -1.  0.  0.]/1.;
bdf2  = [ 3.  -4.  1.  0.]/2.;
bdf3  = [11. -18.  9. -2.]/6.;
ex0   = [0. 1.  0. 0.];
ex1   = [0. 2. -1. 0.];
ex2   = [0. 3. -3. 1.];

u    = zeros(Float64,2);
Ulag = zeros(Float64,2,3);
Rlag = zeros(Float64,2,2);


#rhs  = zeros(Float64,2)       # temporary

# Initial condition
y0 = 0.1
x0 = 0.0

u0 = [x0, y0]
u  = copy(u0)

A  = Matrix{Float64}(I,2,2)
Ainv = copy(A)

α = -0.0          # Negative is unstable
ϵ = 0.2

nsteps = 1000000
t = 0.0
dt = 0.0001

xhist = zeros(Float64,nsteps)

for i in 1:nsteps
  global A,Ainv,t
  global u,Ulag,Rlag
  global xhist

#  global rhs

  t = t + dt;

  bdf = bdf2;
  ext = ex1;
     
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

  if (i<=3) 
    A[1,1] = bdf[1]/dt - ϵ + α
    A[1,2] = 1.0
    A[2,1] = -1.0
    A[2,2] = bdf[1]/dt
  
    Ainv   = inv(A)
  end  

  rhsy = -1.0/3.0*ϵ*(u[1]^3)
  rhsx = 0.0
  rhs1 = [rhsy, rhsx]

  rhs = ext[2]*rhs1 .+ ext[3]*Rlag[1] .+ ext[4]*Rlag[2];
  Rlag[:,2] = Rlag[:,1];
  Rlag[:,1] = rhs1;

  bdlag = 1. /dt*(bdf[2]*u .+ bdf[3]*Ulag[:,1] .+ bdf[4]*Ulag[:,2]); 
  Ulag[:,3] = Ulag[:,2];
  Ulag[:,2] = Ulag[:,1];
  Ulag[:,1] = u;

  rhs = rhs .- bdlag 

  u   = Ainv*rhs
  y   = u[1]
  x   = u[2]

  xhist[i] = x
end


time = range(dt,step=dt,length=nsteps)./π
plot(time,xhist);






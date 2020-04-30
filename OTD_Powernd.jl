#!/usr/bin/julia
# Source code for Exercise 2 of Homework 1


println("OTD comparison with Power method")

import Pkg
#Pkg.add("PyPlot")
#Pkg.add("Blink")

using Blink
using PyPlot
using LinearAlgebra


close("all")


# Temporal discretization
bdf1  = [ 1.  -1.  0.  0.]/1.;
bdf2  = [ 3.  -4.  1.  0.]/2.;
bdf3  = [11. -18.  9. -2.]/6.;
ex0   = [0. 1.  0. 0.];
ex1   = [0. 2. -1. 0.];
ex2   = [0. 3. -3. 1.];


# Create Matrix A
n = 100;
v0 = -0.2 .+ -1. *rand(Float64,n);

v0[1] = -0.01;
v0[2] = -0.05;
v0[3] = -0.1;
v0[4] = -0.15;

n = length(v0);
A = zeros(Float64,n,n);

for i in 1:n
  A[i,i] = v0[i];
end

# Parameters controlling the modal decay rates

dt = 0.01;
Nstep = 20000;
egvupd = 1000;


nmodes = 3;
Vinit  = rand(Float64,n,nmodes);
V      = zeros(Float64,n,nmodes);
Rhs    = zeros(Float64,n,nmodes);
Vlag   = zeros(Float64,n,3,nmodes);
R1lag  = zeros(Float64,n,2,nmodes);
R2lag  = zeros(Float64,n,2,nmodes);

# Power method variables
Vpinit = zeros(Float64,n,nmodes);
Vp     = zeros(Float64,n,nmodes);
Vplag  = zeros(Float64,n,3,nmodes);
Rplag  = zeros(Float64,n,2,nmodes);
Rhsp   = zeros(Float64,n,nmodes);


v = copy(Vinit[:,1]);
#v[1] = 1.e-0;
#v[2] = 1.;
v = v/norm(v);
Vinit[:,1] = v;

# Orthogonalize initial field
for i=2:nmodes
  v = Vinit[:,i];
  alpha = (Vinit[:,1:i-1])'*v;
  w = v - Vinit[:,1:i-1]*alpha;
  w = w/norm(w);
  Vinit[:,i] = w;
end

Vpinit = copy(Vinit);
V      = copy(Vinit);
Vp     = copy(Vpinit);

VNORM = zeros(Float64,Nstep,nmodes);
Evals = zeros(Complex,Nstep,nmodes);
Ermax = zeros(Complex,Nstep);
Epvals = zeros(Complex,Nstep,nmodes);


t     = 0.

h1  = figure(num=1,figsize=[8.,6.]);
   

time = range(dt,step=dt,length=Nstep);

ee0 = eigvals(A);
ee  = sort(ee0,rev=true);
emax1 = ee[1]*ones(Float64,Nstep);
pl1 = plot(time,emax1,linestyle="--")
emax1 = ee[2]*ones(Float64,Nstep);
pl1 = plot(time,emax1,linestyle="--")

ax1 = h1.gca();

for i in 1:Nstep
  global V, Vlag,Rlag,Rhs,Evals,Ermax
  global t, ax2
    
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

  A1     = A;

  Ar = V'*A1*V;
  ee = eigvals(Ar);
  ee2 = sort(ee,rev=true);
  Evals[i,:] = ee2;

#  Ar = Ar + Ar';
  ee = eigvals(Ar);
  Ermax[i] = maximum(ee);

# Clear the previous figure  
#  if (mod(i,egvupd)==0)
#    ax2.cla()
#  end  


# OTD Calculation 
  for j in 1:nmodes
    global V,Vlag,R1lag,R2lag,Rhs
    global pl2,pl3

    global Vp, Vplag,Rplag,Rhsp,Epvals

#   Extrapolate A*x
    rhs     = A1*V[:,j];
    rhs1    = ext[2]*rhs + ext[3]*R1lag[:,1,j] + ext[4]*R1lag[:,2,j];

    R1lag[:,2,j] = R1lag[:,1,j];
    R1lag[:,1,j] = rhs;


    Ax      = A1*V[:,j]; 
    alpha   = V'*Ax; 
    w       = V*alpha;
    rhs2    = (ext[2]*w + ext[3]*R2lag[:,1,j] + ext[4]*R2lag[:,2,j]);

    R2lag[:,2,j] = R2lag[:,1,j];
    R2lag[:,1,j] = w;

    bdlag    = 1. /dt*(bdf[2]*V[:,j] + bdf[3]*Vlag[:,1,j] + bdf[4]*Vlag[:,2,j]); 
    Rhs[:,j] = rhs1 - rhs2 - bdlag;

    Vlag[:,3,j] = Vlag[:,2,j];
    Vlag[:,2,j] = Vlag[:,1,j];
    Vlag[:,1,j] = V[:,j];

    V[:,j]      = copy(Rhs[:,j])*dt/bdf[1];

#   Power Method Calculation
#   Extrapolate A*x
    rhs     = A1*Vp[:,j];
    rhs1    = ext[2]*rhs + ext[3]*Rplag[:,1,j] + ext[4]*Rplag[:,2,j];

    Rplag[:,2,j] = Rplag[:,1,j];
    Rplag[:,1,j] = rhs;

    bdlag   = 1. /dt*(bdf[2]*Vp[:,j] + bdf[3]*Vplag[:,1,j] + bdf[4]*Vplag[:,2,j]); 
    Rhsp    = rhs1 - bdlag;

    Vplag[:,3,j] = Vplag[:,2,j];
    Vplag[:,2,j] = Vplag[:,1,j];
    Vplag[:,1,j] = Vp[:,j];

    Vp[:,j]      = copy(Rhsp)*dt/bdf[1];
    v2           = Vp[:,j];
    Rayleigh     = v2'*A*v2/(v2'*v2);

    Epvals[i,j]  = Rayleigh;

    if (mod(i,egvupd)==0)

      ax1.plot(t,real(Evals[i,j]),marker=".",color="black")
    
      ax1.plot(t,real(Epvals[i,j]),marker=".",color="red")

      ax1.set_xlabel(L"time")
      ax1.set_ylabel(L"\lambda")
      ax1.set_title("Approximated Eigenvalue")
 
      pause(0.0001)

    end  

  end

# Orthogonalize field
# Do we also need to rescale?
  v = V[:,1];
  v = v/norm(v);
  V[:,1] = v;

  for j=2:nmodes
    v = V[:,j];
    alpha = (V[:,1:j-1])'*v;
    w = v - V[:,1:j-1]*alpha;
    w = w/norm(w);
    V[:,j] = w;
  end


# Normalize Vp  
  v = Vp[:,1];
  v = v/norm(v);
  Vp[:,1] = v;

  for j=2:nmodes
    v = Vp[:,j];
    alpha = (Vp[:,1:j-1])'*v;
    w = v - Vp[:,1:j-1]*alpha;
    w = w/norm(w);
    Vp[:,j] = w;
  end

end

h2 = figure(num=2);
for j in 1:nmodes
  err1 = abs.(v0[j] .- Evals[:,j]);
  semilogy(time,err1);
  err2 = abs.(v0[j] .- Epvals[:,j]);
  semilogy(time,err2,linestyle="--");
end  
ax3 = gca();
ax3.set_xlabel(L"time")
ax3.set_ylabel(L"||\epsilon||")

c1 = exp(v0[2]*dt)/exp(v0[1]*dt); 

itr = range(0,step=1,length=Nstep);
conv1  = c1.^itr;
semilogy(time,conv1,linestyle=":");


println("Done.")






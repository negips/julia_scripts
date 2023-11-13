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
v0 = [-0.01 -0.02];

n = length(v0);
A = zeros(Float64,n,n);

for i in 1:n
  A[i,i] = v0[i];
end

# Parameters controlling the modal decay rates

dt = 0.01;
Nstep = 50000;
egvupd = 1000;


nmodes = 1;
Vinit  = rand(Float64,n,nmodes);
V      = zeros(Float64,n,nmodes);
Rhs    = zeros(Float64,n,nmodes);
Vlag   = zeros(Float64,n,3,nmodes);
R1lag  = zeros(Float64,n,2,nmodes);
R2lag  = zeros(Float64,n,2,nmodes);

# Power method variables
Vp     = zeros(Float64,n);
Vplag  = zeros(Float64,n,3);
Rplag  = zeros(Float64,n,2);
Rhsp   = zeros(Float64,n);


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

V   = copy(Vinit);
Vp  = copy(V[:,1]);

VNORM = zeros(Float64,Nstep,nmodes);
Evals = zeros(Complex,Nstep,nmodes);
Ermax = zeros(Complex,Nstep);

Epvals = zeros(Complex,Nstep);


t     = 0.

h1  = figure(num=1,figsize=[18.,6.]);
ax1 = subplot(121);
ax2 = subplot(122);
   

time = range(dt,step=dt,length=Nstep);

ee = eigvals(A);
emax1 = ee[1]*ones(Float64,Nstep);
pl1 = ax1.plot(time,emax1,linestyle="--")
emax1 = ee[2]*ones(Float64,Nstep);
pl1 = ax1.plot(time,emax1,linestyle="--")


for i in 1:Nstep
  global V, Vlag,Rlag,Rhs,Evals,Ermax
  global t, ax2
  global Vp, Vplag,Rplag,Rhsp,Epvals
    
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
  Evals[i,:] = ee;

#  Ar = Ar + Ar';
  ee = eigvals(Ar);
  Ermax[i] = maximum(ee);

# Clear the previous figure  
  if (mod(i,egvupd)==0)
    ax2.cla()
  end  


# OTD Calculation 
  for j in 1:nmodes
    global V,Vlag,R1lag,R2lag,Rhs

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

    vdiff       = V[:,j] .- Vlag[:,j];

    VNORM[i,j]  = norm(V[:,j]);

    if (mod(i,egvupd)==0)
#      ax2.cla()
#      ax2.plot([0., V[1,j]],[0., V[2,j]]);
      pl2 = ax2.arrow(0.,0.,V[1,j],V[2,j],width=0.01,color="black",length_includes_head=true);

      scale = 1.
      pl22 = ax2.arrow(Vlag[1,j],V[2,j],scale*vdiff[1],scale*vdiff[2],width=0.005,color="blue",length_includes_head=true);     
      ax2.set_xlabel(L"x_{1}")
      ax2.set_ylabel(L"x_{2}")
      ax2.set_xlim([0,1.25])
      ax2.set_ylim([0.,1.25])
      ax2.set_title("Approximated Eigenvctor")
 
      
      ax1.plot(t,real(Ermax[i]),marker=".",color="black")
      ax1.set_xlabel(L"time")
      ax1.set_ylabel(L"\lambda")
      ax1.set_title("Approximated Eigenvalue")

    end  

  end

#   Power Method Calculation
#   Extrapolate A*x
    rhs     = A1*Vp;
    rhs1    = ext[2]*rhs + ext[3]*Rplag[:,1] + ext[4]*Rplag[:,2];

    Rplag[:,2] = Rplag[:,1];
    Rplag[:,1] = rhs;

    bdlag    = 1. /dt*(bdf[2]*Vp + bdf[3]*Vplag[:,1] + bdf[4]*Vplag[:,2]); 
    Rhsp     = rhs1 - bdlag;

    Vplag[:,3] = Vplag[:,2];
    Vplag[:,2] = Vplag[:,1];
    Vplag[:,1] = Vp;

    Vp         = copy(Rhsp)*dt/bdf[1];
    Rayleigh   = Vp'*A*Vp/(Vp'*Vp);

    Epvals[i]  = Rayleigh;

    if (mod(i,egvupd)==0)
      pl3 = ax2.arrow(0.,0.,Vp[1],Vp[2],width=0.01,color="red",length_includes_head=true);
      ax1.plot(t,real(Rayleigh),marker=".",color="red")

      pause(0.001)

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
  Vp = Vp/norm(Vp);

end

h2 = figure(num=2);
err1 = abs.(v0[1] .- Evals);
semilogy(time,err1);
err2 = abs.(v0[1] .- Epvals);
semilogy(time,err2);
ax3 = gca();
ax3.set_xlabel(L"time")
ax3.set_ylabel(L"||\epsilon||")



println("Done.")






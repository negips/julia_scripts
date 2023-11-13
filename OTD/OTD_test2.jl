#!/usr/bin/julia
# Source code for Exercise 2 of Homework 1

println("Modified OTD test in a Jordan form matrix")

import Pkg
#Pkg.add("PyPlot")
#Pkg.add("Blink")

using Blink
using PyPlot
using LinearAlgebra


# Temporal discretization
bdf1  = [ 1.  -1.  0.  0.]/1.;
bdf2  = [ 3.  -4.  1.  0.]/2.;
bdf3  = [11. -18.  9. -2.]/6.;
ex0   = [0. 1.  0. 0.];
ex1   = [0. 2. -1. 0.];
ex2   = [0. 3. -3. 1.];


# Create Matrix A
v = [-0.5, -0.5, -0.2, -0.2, -0.1, -0.1];

n = length(v);
A = zeros(Float64,n,n);

for i in 1:n
  A[i,i] = v[i];
end

# Parameters controlling the various non-modal growths.
# Assuming decay rates are constant.
nm10 = 20.0e-0;
nm1a = 0.;

nm20 = 5.0e-0;
nm2a = 0.;

nm30 = 3.0e-0;
nm3a = 0.0;

# Parameters controlling the modal decay rates
w10 = -1.0;
w1a = 0.;

w20 = -0.1;
w2a = 0.;

w30 = -0.01;
w3a = 0.0;

omega = 0.1;
dt = 0.01;
Nstep = 50000;


nmodes = 4;
Vinit  = rand(Float64,n,nmodes);
V      = zeros(Float64,n,nmodes);
Vlag   = zeros(Float64,n,3,nmodes);
R1lag  = zeros(Float64,n,2,nmodes);
R2lag  = zeros(Float64,n,2,nmodes);

v = 0.0001*Vinit[:,1];
#v[1] = 1.0;
#v[2] = 1.2;
v = v/norm(v);
Vinit[:,1] = v;

# v = 0.0001*Vinit[:,2];
# v[1] = -1.0;
# v[2] = 3.7;
# v = v/norm(v);
# Vinit[:,2] = v;

# Orthogonalize initial field
for i=2:nmodes
  v = Vinit[:,i];
  alpha = (Vinit[:,1:i-1])'*v;
  w = v - Vinit[:,1:i-1]*alpha;
  w = w/norm(w);
  Vinit[:,i] = w;
end

V   = copy(Vinit);
Rhs = 0*copy(V);

VNORM = zeros(Float64,Nstep,nmodes);
Evals = zeros(Complex,Nstep,nmodes);
Ermax = zeros(Complex,Nstep);

t     = 0.

for i in 1:Nstep
  global V, Vlag,Rlag,Rhs,Evals,Ermax
  global t
     
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
  
  A1[1,1] = w10  + w1a*sin(omega*t);
  A1[2,2] = w10  + w1a*sin(omega*t);
  A1[1,2] = nm10 + nm1a*sin(omega*t);

  A1[3,3] = w20  + w2a*sin(omega*t);
  A1[4,4] = w20  + w2a*sin(omega*t);
  A1[3,4] = nm20 + nm2a*sin(omega*t);

  A1[5,5] = w30  + w3a*sin(omega*t);
  A1[6,6] = w30  + w3a*sin(omega*t);
  A1[5,6] = nm30 + nm3a*sin(omega*t);


  Ar = V'*A1*V;
  ee = eigvals(Ar);
  Evals[i,:] = ee;

  Ar = Ar + Ar';
  ee = eigvals(Ar);
  Ermax[i] = maximum(ee);

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
    rhs2    = 1*(ext[2]*w + ext[3]*R2lag[:,1,j] + ext[4]*R2lag[:,2,j]);


    R2lag[:,2,j] = R2lag[:,1,j];
    R2lag[:,1,j] = w;

    bdlag    = 1. /dt*(bdf[2]*V[:,j] + bdf[3]*Vlag[:,1,j] + bdf[4]*Vlag[:,2,j]); 
    Rhs[:,j] = rhs1 - rhs2 - bdlag;

    Vlag[:,3,j] = Vlag[:,2,j];
    Vlag[:,2,j] = Vlag[:,1,j];
    Vlag[:,1,j] = V[:,j];

    V[:,j]      = copy(Rhs[:,j])*dt/bdf[1];

    VNORM[i,j]  = norm(V[:,j]);

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


end

#deleteat!(err_hist,1)

# Not sure what this is called. Loading a module?
# using Plotly

#j = [i for i in 1:Niters];
#
time = range(dt,step=dt,length=Nstep);
layout = ["yaxis" => ["type"=>"log", "autorange" => true] ];
# pl = plot(err_hist, ["layout" => layout])
pl = plot(time,real.(Evals))
pl2 = plot(time,real.(Ermax),linestyle=":")
#ylabel("|exp(At)|");
#xlabel("time");


#AA = A + A';
#ee = eigvals(AA);

#emax = maximum(ee).*ones(Float64,Nstep);
#pl3 = plot(time,emax,linestyle="--")

emax1 = maximum(ee).*ones(Float64,Nstep);
emax2 = maximum(ee).*ones(Float64,Nstep);
emax3 = maximum(ee).*ones(Float64,Nstep);

for i in 1:Nstep
  global emax1,emax2,emax3

  t = time[i];

  A1 = A[1:2,1:2];
  A1[1,1] = w10  + w1a*sin(omega*t);
  A1[2,2] = w10  + w1a*sin(omega*t);
  A1[1,2] = nm10 + nm1a*sin(omega*t);

  ee = eigvals(A1 + A1');
  emax1[i] = maximum(ee);

  A2 = A[3:4,3:4];
  A2[1,1] = w20  + w2a*sin(omega*t);
  A2[2,2] = w20  + w2a*sin(omega*t);
  A2[1,2] = nm20 + nm2a*sin(omega*t);

  ee = eigvals(A2 + A2');
  emax2[i] = maximum(ee);
  
  A3 = A[5:6,5:6];
  A3[1,1] = w30  + w3a*sin(omega*t);
  A3[2,2] = w30  + w3a*sin(omega*t);
  A3[1,2] = nm30 + nm3a*sin(omega*t);

  ee = eigvals(A3 + A3');
  emax3[i] = maximum(ee);

end

pl4 = plot(time,emax1,linestyle="--")

pl5 = plot(time,emax2,linestyle="--")

pl6 = plot(time,emax3,linestyle="--")


println("Done.")






#!/usr/bin/julia
# Source code for Exercise 2 of Homework 1

println("OTD test in a Jordan form matrix")

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
v = [-1.0, -1.0, -0.1, -0.1, -0.01, -0.01];

n = length(v);
A = zeros(Float64,n,n);

for i in 1:n
  A[i,i] = v[i];
end

epsilon = 1.0e-1;
A[1,2] = epsilon;

epsilon = 1.00e-0;
A[3,4] = epsilon;

epsilon = 0.1e-0;
A[5,6] = epsilon;


nmodes = 2;
Vinit  = rand(Float64, (n,nmodes));
V      = zeros(Float64, (n,nmodes));
Vlag   = zeros(Float64,n,3,nmodes);
Rlag   = zeros(Float64,n,2,nmodes);

v = Vinit[:,1];
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

V   = Vinit;
Rhs = 0*copy(V);

dt = 0.01;
Nstep = 50000;

VNORM = zeros(Float64,Nstep,nmodes);
Evals = zeros(Complex,Nstep,nmodes);
Ermax = zeros(Complex,Nstep);

for i in 1:Nstep
  global V, Vlag,Rlag,Rhs,Evals,Ermax
     
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

  Ar = V'*A*V;
  ee = eigvals(Ar);
  Evals[i,:] = ee;

  Ar = Ar + Ar';
  ee = eigvals(Ar);
  Ermax[i] = maximum(ee);

  for j in 1:nmodes
    global V, Vlag,Rlag,Rhs

#   Extrapolate A*x
    rhs1    = A*(ext[2]*V[:,j] + ext[3]*Vlag[:,1,j] + ext[4]*Vlag[:,2,j]);

    Ax      = A*V[:,j]; 
    alpha   = V'*Ax; 
    w       = V*alpha;
    rhs2    = 1*(ext[2]*w + ext[3]*Rlag[:,1,j] + ext[4]*Rlag[:,2,j]);


    Rlag[:,2,j] = Rlag[:,1,j];
    Rlag[:,1,j] = w;

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


AA = A + A';
ee = eigvals(AA);

emax = maximum(ee).*ones(Float64,Nstep);
#pl3 = plot(time,emax,linestyle="--")


A1 = A[1:2,1:2] + A[1:2,1:2]';
ee1 = eigvals(A1);

emax = maximum(ee1).*ones(Float64,Nstep);
pl4 = plot(time,emax,linestyle="--")

A2 = A[3:4,3:4] + A[3:4,3:4]';
ee2 = eigvals(A2);

emax = maximum(ee2).*ones(Float64,Nstep);
pl5 = plot(time,emax,linestyle="--")

A3 = A[5:6,5:6] + A[5:6,5:6]';
ee3 = eigvals(A3);

emax = maximum(ee3).*ones(Float64,Nstep);
pl6 = plot(time,emax,linestyle="--")


println("Done.")






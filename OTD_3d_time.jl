#!/usr/bin/julia
# Source code for Exercise 2 of Homework 1

println("Visualizing the changing approximated eigenvectors/eigenvalues")
println("from the OTD method in a 3x3 system with 2 OTD modes")

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

close("all")

# Create Matrix A
v0 = [-0.001 -0.020 -0.015];
v3a    = 0.004;
Omega  = 0.002;

n = length(v0);
A = zeros(Float64,n,n);

for i in 1:n
  A[i,i] = v0[i];
  for j in i+1:n
    tmp = rand(Float64)
    A[i,j] = tmp;
  end  
end


# Parameters controlling the modal decay rates

dt = 0.05;
Nstep = 150000;
egvupd = 1000;


nmodes = 2;
Vinit  = rand(Float64,n,nmodes);
V      = zeros(Float64,n,nmodes);
Vlag   = zeros(Float64,n,3,nmodes);
R1lag  = zeros(Float64,n,2,nmodes);
R2lag  = zeros(Float64,n,2,nmodes);

v = copy(Vinit[:,1]);
#v[1] = 1.e-3;
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
Rhs = 0*copy(V);

VNORM = zeros(Float64,Nstep,nmodes);
Evals = zeros(Complex,Nstep,nmodes);
Ermax = zeros(Complex,Nstep);

t     = 0.

v1    = [0, 1.];
v2    = [1., 0.];

h1  = figure(num=1,figsize=[18.,6.]);
ax1 = subplot(121);
ax2 = subplot(122,projection="3d");

   

time = range(dt,step=dt,length=Nstep);

ee0 = eigvals(A);
ee  = sort(ee0,rev=true);

emax1 = zeros(Float64,Nstep);
emax2 = zeros(Float64,Nstep);
emax3 = zeros(Float64,Nstep);

for i in 1:Nstep
  global emax1,emax2,emax3
  emax1[i] = v0[1];
  emax2[i] = v0[2];
  emax3[i] = v0[3] + v3a*sin(Omega*time[i]);
end  

pl1 = ax1.plot(time,emax1,linestyle="--")
pl2 = ax1.plot(time,emax2,linestyle="--")
pl3 = ax1.plot(time,emax3,linestyle="--")

for i in 1:Nstep
  global V, Vlag,Rlag,Rhs,Evals,Ermax
  global t, ax2, A1
     
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

  A1     = copy(A);
  A1[1,1] = v0[1];
  A1[2,2] = v0[2];
  A1[3,3] = v0[3] + v3a*sin(Omega*time[i]);

  Ar = V'*A1*V;
  ee = eigvals(Ar);
  Evals[i,:] = ee;

#  Ar = Ar + Ar';
  ee = eigvals(Ar);
  for j=1:nmodes
    Evals[i,j] = ee[j];
  end
  ee = eigvals(Ar + Ar');
  Ermax[i]   = maximum(ee);

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

    if (mod(i,egvupd)==0)
      if (j==1)
        ax2.cla()
      end

      ax2.plot([0., V[1,j]], [0., V[2,j]], [0., V[3,j]]);
#      pl2 = ax2.arrow(0.,0.,0.,V[1,j],V[2,j],V[3,j],width=0.02,color="black",length_includes_head=true);
      ax2.set_xlabel(L"x_{1}")
      ax2.set_ylabel(L"x_{2}")
      ax2.set_zlabel(L"x_{3}")
     
      ax2.set_xlim([-1.25,1.25])
      ax2.set_ylim([-1.25,1.25])
      ax2.set_zlim([-1.25,1.25])
     
      ax2.set_title("Approximated Basis")
     
      if (j==1)
#        ax1.plot(t,real(Ermax[i]),marker=".",color="gray")
      end  

      ax1.plot(t,real(Evals[i,j]),marker=".",color="black")
      ax1.set_xlabel(L"time")
      ax1.set_ylabel(L"\lambda")
      ax1.set_title("Approximated Eigenvalues")

      if j==2
        pause(0.00001)
      end  
    end  
#    h1.canvas.draw() # Update the figure

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


println("Done.")






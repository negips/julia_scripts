#!/usr/bin/julia
# Source code for Exercise 2 of Homework 1


println("Visualizing the changing approximated eigenvector/eigenvalue")
println("from the OTD method in a 2x2 system, for a single OTD mode")

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

close("all")

dt     = 0.001;
nsteps = 200000;
figupd = 100;

a1    = 2.;
a2    = 2.;
ϵ     = 0.05;
b     = 20.;

time = range(dt,step=dt,length=nsteps);

# Random.seed!(17);           # 22, 30, 31, 37 

Xinit = rand(Float64,n);
X     = zeros(Float64,n);
Xlag  = zeros(Float64,n,3);
Rlag  = zeros(Float64,n,2);

X     = copy(Xinit);
#X     = [2.0; 1.; 0.1]
X[3]  = -1.;


Rhs   = 0*copy(X);

t     = 0.

h1  = figure(num=1,figsize=[8.,6.]);
#ax  = PyPlot.axes3D();
#ax2 = subplot(122,projection="3d");

a1d   = ones(Float64,1);
x1h   = zeros(Float64,nsteps);
x2h   = zeros(Float64,nsteps);
x3h   = zeros(Float64,nsteps);


for i in 1:nsteps
  global X, Xlag,Rlag,Rhs
  global t, a1d
  global x1h,x2h,x3h
  global pl,sc

     
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

# .  
# x1 = -a1x1 + ϵx2 + bx3
#
# .
# x2 = (1/ϵ)x1 - a2x2
#
# .
# x3 = bx3*(1/√(x1^2 + x2^2))


# Extrapolate A*x
#
  x1        = X[1];
  x2        = X[2];
  x3        = X[3];
 
  x1r       = -a1*x1 + ϵ*x2 + b*x3;
  rhs1      = ext[2]*x1r + ext[3]*Rlag[1,1] + ext[4]*Rlag[1,2];

  x2r       = (1/ϵ)*x1 - a2*x2;
  rhs2      = ext[2]*x2r + ext[3]*Rlag[2,1] + ext[4]*Rlag[2,2];

  x3r       = b*(1/sqrt(x1^2 + x2^2) - 1.)*x3;
  rhs3      = ext[2]*x3r + ext[3]*Rlag[3,1] + ext[4]*Rlag[3,2];
 
  rhs       = [x1r; x2r; x3r];
  rhsall    = [rhs1; rhs2; rhs3];

  Rlag[:,2] = Rlag[:,1];
  Rlag[:,1] = rhs;

  bdlag     = 1. /dt*(bdf[2]*X + bdf[3]*Xlag[:,1] + bdf[4]*Xlag[:,2]); 
  Rhs       = rhsall - bdlag;

  Xlag[:,3] = Xlag[:,2];
  Xlag[:,2] = Xlag[:,1];
  Xlag[:,1] = X;

  X         = Rhs*dt/bdf[1];

  xdiff     = X - Xlag[:,1];
  
  x1p       = X[1]*a1d;
  x2p       = X[2]*a1d;
  x3p       = X[3]*a1d;
#

  x1h[i]     = X[1];
  x2h[i]     = X[2];
  x3h[i]     = X[3];
 

  if (mod(i,figupd)==0)
    if i>figupd
      pl[1].remove();
      sc.remove();
    end

    nhis = 8000;
    if i<nhis
      j=1:i;
    else
      j=(i-nhis+1):i;
    end

    pl = plot3D(x1h[j],x2h[j],x3h[j],color=rgba0);
    sc = scatter3D(x1p,x2p,x3p,color="red")
    ax = gca();
    ax.set_xlabel(L"x_{1}")
    ax.set_ylabel(L"x_{2}")
    ax.set_zlabel(L"x_{3}")

#    ax.set_xlim([-2,2])
#    ax.set_ylim([-2,10])
#    ax.set_zlim([-1.,1.])

    pause(0.001)
  end    
        


end

println("Done.")






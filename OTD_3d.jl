#!/usr/bin/julia
# Source code for Exercise 2 of Homework 1

println("Visualizing the changing approximated eigenvectors/eigenvalues")
println("from the OTD method in a 3x3 system with 2 OTD modes")

import Pkg
#Pkg.add("PyPlot")
#Pkg.add("Blink")

using Blink
using PyPlot
using LinearAlgebra,Random


nonnormal = true;
ifoptimal = false;


# Temporal discretization
bdf1  = [ 1.  -1.  0.  0.]/1.;
bdf2  = [ 3.  -4.  1.  0.]/2.;
bdf3  = [11. -18.  9. -2.]/6.;
ex0   = [0. 1.  0. 0.];
ex1   = [0. 2. -1. 0.];
ex2   = [0. 3. -3. 1.];

close("all")

# Create Matrix A
v0 = [-0.005 -0.020 -0.015];

n = length(v0);
A = zeros(Float64,n,n);

for i in 1:n
  A[i,i] = v0[i];
end

if nonnormal
  global A
#  A[1,2] = -0.05;
 
  A[2,3] = 0.06;
end  


# Parameters controlling the modal decay rates

dt = 0.05;
Nstep = 20000;
egvupd = 200;


Evals = zeros(Complex,Nstep,nmodes);
Ermax = zeros(Complex,Nstep);


h1  = figure(num=1,figsize=[18.,6.]);
ax1 = subplot(121);
ax2 = subplot(122,projection="3d");
ax2.set_xlabel(L"x_{1}")
ax2.set_ylabel(L"x_{2}")
ax2.set_zlabel(L"x_{3}")

ax2.set_xlim([-1.25,1.25])
ax2.set_ylim([-1.25,1.25])
ax2.set_zlim([-1.25,1.25])

ax2.set_title("Approximated Basis")

cm    = get_cmap("tab10");
rgba0 = cm(0); 
rgba1 = cm(1); 
rgba2 = cm(2); 
rgba3 = cm(3); 


time = range(dt,step=dt,length=Nstep);

F0  = eigen(A);
ee  = sortperm(F0.values,rev=true);

emax1 = F0.values[ee[1]].*ones(Float64,Nstep);
emax2 = F0.values[ee[2]].*ones(Float64,Nstep);
emax3 = F0.values[ee[3]].*ones(Float64,Nstep);

Fe    = eigen(0.5*(A + A'));
ee2   = sortperm(Fe.values,rev=true);
emax4 = Fe.values[ee2[1]].*ones(Float64,Nstep);
emax5 = Fe.values[ee2[2]].*ones(Float64,Nstep);

ple1 = ax1.plot(time,emax1,linestyle="--",color=rgba0)
ple2 = ax1.plot(time,emax2,linestyle="--",color=rgba1)
ple3 = ax1.plot(time,emax3,linestyle="--",color=rgba2)
ple4 = ax1.plot(time,emax4,linestyle=":",color="gray")
ple5 = ax1.plot(time,emax5,linestyle=":",color="gray")

ax1.set_xlabel(L"time")
ax1.set_ylabel(L"\lambda")
ax1.set_title("Approximated Eigenvalues")


Random.seed!(10);       # 10

nmodes = 2;
Vinit  = rand(Float64,n,nmodes);
V      = zeros(Float64,n,nmodes);
Vlag   = zeros(Float64,n,3,nmodes);
R1lag  = zeros(Float64,n,2,nmodes);
R2lag  = zeros(Float64,n,2,nmodes);

#Vinit[:,1] = 0.1*F0.vectors[:,ee[2]] +  1.0*F0.vectors[:,ee[1]]
#Vinit[:,2] = 0.5*F0.vectors[:,ee[1]] +  0.0*F0.vectors[:,ee[3]]

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


t     = 0.

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

  if (ifoptimal)
    A1   = 0.5*(A + A');
  else  
    A1     = copy(A);
  end  
#  A1[1,1] = v0[1];
#  A1[2,2] = v0[2];
#  A1[3,3] = v0[3] + v3a*sin(Omega*time[i]);

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
    global pl5,pl6,pl7,pl8

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

    if (mod(i,egvupd)==0)
      if (j==1)
#        ax2.cla()
      end

      if i>egvupd
        if j==1
          pl5[1].remove()
          pl6[1].remove()
        else
          pl7[1].remove()
          pl8[1].remove()
        end  
      end


#      pl2 = ax2.arrow(0.,0.,0.,V[1,j],V[2,j],V[3,j],width=0.02,color="black",length_includes_head=true);
#      ax2.set_xlabel(L"x_{1}")
#      ax2.set_ylabel(L"x_{2}")
#      ax2.set_zlabel(L"x_{3}")
#     
#      ax2.set_xlim([-1.25,1.25])
#      ax2.set_ylim([-1.25,1.25])
#      ax2.set_zlim([-1.25,1.25])
#     
#      ax2.set_title("Approximated Basis")
     
      if (j==1)
#        ax1.plot(t,real(Ermax[i]),marker=".",color="gray")
      end  

      if j==1
        pl5 = ax2.plot([0., V[1,j]], [0., V[2,j]], [0., V[3,j]],color=rgba0);
        pl6 = ax1.plot(time[1:i],real(Evals[1:i,j]),color="black");
      else
        pl7 = ax2.plot([0., V[1,j]], [0., V[2,j]], [0., V[3,j]],color=rgba1);
        pl8 = ax1.plot(time[1:i],real(Evals[1:i,j]),color="black");
      end  

      if j==nmodes
        pause(0.00001)
      end  
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


end


println("Done.")






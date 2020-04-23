#!/usr/bin/julia
# Source code for Homework 3

println("Homework 2, Exercise 6")

#import Pkg
#Pkg.add("PyPlot")
#Pkg.add("PyCall")
#Pkg.add("Blink")

using Blink,Random
using PyPlot,PyCall,LinearAlgebra,SparseArrays
using LaTeXStrings,TimerOutputs,MAT

include("mycg.jl")
include("mygmres.jl")
include("GS_npass.jl")

close("all")      # close all open figures

const TO = TimerOutput();

ifsaveerr=true;
ifsavetime=false;
lgfs = 8;
lafs = 14;
ifmod = false;    # if we use modified GS
ngs = 1;          # No of GS passes

cols = ["b","r","k","g"]

vars = matread("Bwedge2.mat");
A = Matrix(vars["B"]);
b = Complex.(vars["b"]);

r,c = size(A);

xexact = A\b;

tol = 1.0e-15;

range = [-13]; #collect(-2:-1:-2)
tol_range = 10. .^range;
ntol = length(tol_range);

# Dummy for timing
ifcomplex=true;
ifcgn = true;
maxits = 1;
x1,X1,P1,R1,residuals1,xerr1 = @timeit TO "MyCGN.dummy" mycg(A,b,tol,maxits,xexact,ifcgn,ifcomplex);


# CGN
ifcgn = true;
maxits = 120;
cg_time = zeros(Float64,ntol,1);

for i=1:ntol
  global rc    

  tol = tol_range[i];
  xc,Xc,Pc,Rc,rc,xerrc = @timeit TO "MyCGN.$i" mycg(A,b,tol,maxits,xexact,ifcgn,ifcomplex);

  cg_time[i]=TimerOutputs.time(TO["MyCGN.$i"])*10^-9;
  
  if (ifcgn)
    for i=1:length(rc)
      rc[i]=norm(A*Xc[:,i]-b);
    end
  end
end
println("CNG done.")
  

# Dummy for timing
ngs=1;
maxits=1;
ifcomplex=true;
x1,Q1,H1,r1,xerr1 = @timeit TO "MyGMRES.dummy" mygmres(A,b,ngs,ifmod,tol,maxits,xexact,ifcomplex);


# GMRES
gm_time = zeros(Float64,ntol,1);
for i=1:ntol
  tol = tol_range[i];

  global rg

  ngs=1;
  ifmod=false;
  maxits=120;
  xg,Qg,Hg,rg,xerrg = @timeit TO "MyGMRES.$i" mygmres(A,b,ngs,ifmod,tol,maxits,xexact,ifcomplex);
  gm_time[i]=TimerOutputs.time(TO["MyGMRES.$i"])*10^-9;
end  

# Plotting

# Error plot

if (ifsaveerr)
  niters = length(rg);
  iters = collect(1:niters);
  str3=latexstring("\$||Ax - b||\$")
  str1="GMRES";
  pgm = semilogy(iters,rg,color=cols[1],linestyle="-",label=str1);
  
  str2="CGN"
  niters = length(rc);
  iters  = collect(1:niters);
  pcg = semilogy(iters,rc,color=cols[2],linestyle="-",label=str2);
  ylabel(str3,fontsize=lafs)
  xlabel("Iterations",fontsize=lafs)
  legend(fontsize=lgfs);

  fname="ex6_err.png"
#  savefig(fname,dpi=150)
end


## Convergence
evals = eigvals(A);

er = real(evals);
ei = imag(evals);

ind = sortperm(er,by=abs);
ers = er[ind];
eis = ei[ind];
neig = length(ers);

rmax = ers[neig];
rmin = ers[1];          # Skipping the smallest one

cen = (rmax+rmin)/2;
rad = 1.0*(rmax-cen);

r2 = 1.0*(eis[niters]-eis[1])/2;
if r2>rad
  rad=r2;
end  

cfac1 = rad/cen;
theta = collect(0.:0.00001:2*pi);
x1=cen .+ rad*cos.(theta);
y1=0. .+ rad*sin.(theta);

niters = length(rg);
ik = collect(1:niters+3);
conv1 = 1. *rg[1]*(cfac1).^ik;

# Deflated
rmax = ers[neig];
rmin = ers[2];          # Skipping the smallest one
cen = (rmax+rmin)/2;
rad = 1.05*(rmax-cen);
r2 = 1.05*(eis[niters]-eis[1])/2;
if r2>rad
  rad=r2;
end  
cfac2 = rad/cen;
x2=cen .+ rad*cos.(theta);
y2=0. .+ rad*sin.(theta);

conv2 = 500. *rg[1]*(cfac2).^ik;

if ifsaveerr
  pl1 = semilogy(conv1,color=cols[1],linestyle=":",label="Theoretical GMRES")
  pl2 = semilogy(conv2,color=cols[1],linestyle="-.",label="Theoretical GMRES Deflated")
end  

AA = A'*A;
evals2 = eigvals(AA);
er2 = real(evals2);

ind = sortperm(er2,by=abs);
ers2 = er2[ind];

k1 = ers2[neig]/ers2[1];
k2 = ers2[neig]/ers2[2];

cfac1 = (sqrt(k1)-1)/(sqrt(k1)+1);
cfac2 = (sqrt(k2)-1)/(sqrt(k2)+1);

conv1 = 1. *rc[1]*(cfac1).^ik;
conv2 = 1000. *rc[1]*(cfac2).^ik;

if ifsaveerr
  pl3 = semilogy(conv1,color=cols[2],linestyle=":",label="Theoretical CG")
  pl4 = semilogy(conv2,color=cols[2],linestyle="-.",label="Theoretical CG Deflated")

  legend(fontsize=lgfs);
  fname="ex6_err.png"
  savefig(fname,dpi=150)
end  


## Eigenvalues
figure()
fm2 = plt.get_current_fig_manager();
fm2.window.wm_geometry("800x600")
str4=latexstring("\$\\Lambda(A)\$")       
plot(er,ei,linestyle="none",marker=".",markersize=6,label=str4)
plot(x1,y1,linewidth=1.5,color="k",linestyle="-")
plot(x2,y2,linewidth=1.5,color="k",linestyle="--")
ax=gca();
#ax[:set_position]([0.15, 0.15, 0.8, 0.7])
ax.set_position([0.15, 0.15, 0.8, 0.7])
ylabel(L"\lambda_{i}",fontsize=lafs)
xlabel(L"\lambda_{r}",fontsize=lafs)
tight_layout
legend(fontsize=lgfs)
fname="ex6_spectraA.png"
savefig(fname,dpi=150)

# Eigenvalues of A'*A
figure()
fm2 = plt.get_current_fig_manager();
fm2.window.wm_geometry("800x600")
str4=latexstring("\$\\Lambda(A'A)\$")       
plot(ers2,linestyle="none",marker=".",markersize=6,label=str4)
ax=gca();
#ax[:set_position]([0.15, 0.15, 0.8, 0.7])
ax.set_position([0.15, 0.15, 0.8, 0.7])
ylabel(L"\lambda_{r}",fontsize=lafs)
tight_layout
legend(fontsize=lgfs)
fname="ex6_spectraAA.png"
savefig(fname,dpi=150)

# Timing plot
if (ifsavetime)
  figure()
  pgm2 = semilogy(gm_time,tol_range,color=cols[1],linestyle="--",marker=".",label=str1);
  pcg2 = semilogy(cg_time,tol_range,color=cols[2],linestyle="--",marker=".",label=str2);
  ylabel(str3,fontsize=lafs)
  xlabel("Time(s)",fontsize=lafs)
  legend(fontsize=lgfs);

  fname="ex6_timing.png"
#  savefig(fname,dpi=150)
end


show(TO, allocations = false, compact = true)




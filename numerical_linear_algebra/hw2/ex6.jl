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
ifsavetime=true;
ifmod = false;    # if we use modified GS
ngs = 1;          # No of GS passes

cols = ["b","r","k","g"]

vars = matread("Bwedge2.mat");
A = Matrix(vars["B"]);
b = Complex.(vars["b"]);

r,c = size(A);

xexact = A\b;

tol = 1.0e-15;

range = collect(-2:-1:-13)
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
niters = length(rg);
iters = collect(1:niters);
str3=latexstring("\$||Ax - b||\$")
str1="GMRES";
pgm = semilogy(iters,rg,color=cols[1],linestyle="--",marker=".",label=str1);

str2="CGN"
niters = length(rc);
iters  = collect(1:niters);
pcg = semilogy(iters,rc,color=cols[2],linestyle="--",marker=".",label=str2);
ylabel(str3,fontsize="14")
xlabel("Iterations",fontsize="14")
legend(fontsize=14);

if (ifsaveerr)
  fname="ex6_err.png"
  savefig(fname,dpi=150)
end


# Timing plot
figure()
pgm2 = semilogy(gm_time,tol_range,color=cols[1],linestyle="--",marker=".",label=str1);
pcg2 = semilogy(cg_time,tol_range,color=cols[2],linestyle="--",marker=".",label=str2);
ylabel(str3,fontsize="14")
xlabel("Time(s)",fontsize="14")
legend(fontsize=14);

if (ifsavetime)
  fname="ex6_timing.png"
  savefig(fname,dpi=150)
end

show(TO, allocations = false, compact = true)




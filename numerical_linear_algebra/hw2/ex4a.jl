#!/usr/bin/julia
# Source code for Homework 3

println("Homework 2, Exercise 4a")

#import Pkg
#Pkg.add("PyPlot")
#Pkg.add("PyCall")
#Pkg.add("Blink")

using Blink,Random
using PyPlot,PyCall,LinearAlgebra,SparseArrays
using LaTeXStrings,TimerOutputs

include("mycg.jl")

const TO = TimerOutput();

close("all")      # close all open figures

ifmod = false;    # if we use modified GS
ngs = 1;          # No of GS passes

errtol = 1.0e-10; # error tolerance for sub diagonal elements

cols = ["b","r","k","g"];
lafs=20;
lgfs=16;
ifploterr = true;

ud = Vector{Float64}(undef,7);
fill!(ud,1.0);
ud[2]=-1;
md = Vector{Float64}(undef,8);
fill!(md,2.0);

A = SymTridiagonal(md, ud);

b = 0. .*md;
b[1] = 1.;
b[2] = 1.;
 
xexact = @timeit TO "Backslash" A\b;
resid = norm(A*xexact - b);
println("Exact resid=$resid")

tol = 1.0e-14;
maxits = 4;
# Dummy so I get timings right  
x1,X1,P1,R1,r1,xerr1 = mycg(A,b,tol,maxits,xexact); 
x,X,P,R,r,xerr = @timeit TO "MyCG" mycg(A,b,tol,maxits,xexact);

resid=r[end];
println("CG resid=$resid for m=$maxits iterations")

str1=latexstring("\$||x - x_{*}||;\$")
str2=latexstring("\$||Ax - b||;\$")

if (ifploterr)
  pl1 = semilogy(xerr,color=cols[1],label=str1);
  pl2 = semilogy(r,color=cols[1],linestyle="--", label=str2);
  ylabel(L"\varepsilon",fontsize=lafs)
  xlabel("Iterations",fontsize=lafs)
  fm1 = plt.get_current_fig_manager();
  fm1.window.wm_geometry("1200x800")
end  

if (ifploterr)
  legend(loc="best",fontsize=lgfs)
  fname="ex4a_err.png"
  savefig(fname)
end

TO

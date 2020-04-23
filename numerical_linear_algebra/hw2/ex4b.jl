#!/usr/bin/julia
# Source code for Homework 3

println("Homework 2, Exercise 4b")

#import Pkg
#Pkg.add("PyPlot")
#Pkg.add("PyCall")
#Pkg.add("Blink")

using Blink,Random
using PyPlot,PyCall,LinearAlgebra,SparseArrays
using LaTeXStrings

include("mycg.jl")

close("all")      # close all open figures

ifmod = false;    # if we use modified GS
ngs = 1;          # No of GS passes

cols = ["b","r","k","g"]

m = 8;
maxits = 4;
md = 2. *ones(m);
ud = 1. *ones(m-1);
ud[2] = -1.;
A = SymTridiagonal(md,ud)
b = zeros(8);
b[1] = 1.;
b[2] = 1.;
# b[end] = 1.;


# alpha = 10;
# Random.seed!(50);
# A = rand(Float64,m,m);
# A = A + A' + alpha*Matrix(1.0I,m,m);      # Symmetric
# A = A'*A;                                 # SPD
# A = A/norm(A);
# b = rand(Float64,m);
# 
xexact = A\b;
# 

invA = inv(Matrix(A));

tol = 1.0e-14;
# #x,Q,H,r,xerr = mygmres(A,b,ngs,ifmod,tol,maxits,xexact);
x,r,p,xerr,res,R,X,P = mycg(A,b,tol,maxits,xexact);

niters = length(res);
iters = collect(1:niters);
str2=latexstring("\$||Ax - b||_{A^{-1}}\$")
# 
# f2 = figure();
# pl1 = semilogy(xerr,color=cols[1],label=str1);
pl2 = semilogy(iters,res,color=cols[1],linestyle="--",marker=".",label=str2);
ylabel(L"\varepsilon",fontsize="14")
xlabel("Iterations",fontsize="14")
legend();

fm = plt[:get_current_fig_manager]();
fm[:window][:wm_geometry]("800x600")
ax=gca();
#ax[:set_position]([0.15, 0.15, 0.8, 0.7])

res[end]

# 
# #plt[:figure](1)
# fm1 = plt[:get_current_fig_manager]();
# fm1[:window][:wm_geometry]("800x600")



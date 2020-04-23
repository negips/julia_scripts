using LinearAlgebra, Random
using Plots, MatrixDepot, SparseArrays, QuadGK

include("specificGramian.jl")

A = matrixdepot("neumann",40); A=A-A'; A=A/(2*norm(A,1));
B = sprand(Float64,size(A,1),size(A,1),0.05);
objTol = 1e-16;
tau = 1;

P_sp     = @time specificGramian(A,B,tau,objTol);

Pfun(t)  = exp(t*Matrix(A'))*B*exp(t*Matrix(A));
P_naive,errI  = @time quadgk(Pfun,0,tau);

errTrun = norm(P_sp-P_naive);

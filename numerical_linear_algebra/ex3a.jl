
#using Pkg
#Pkg.add("MatrixDepot")
using MatrixDepot,Random,TimerOutputs,PyPlot
using LinearAlgebra # for: eig, norm, etc

const to = TimerOutput();

include("arnoldi.jl")
include("GS_npass.jl")
include("GS_modified.jl")

println("Homework 1, Exercise 3a")

# Create Matrix A
nn=10;
Random.seed!(0)
A=matrixdepot("wathen",nn,nn);
rw,cl = size(A);

b=randn(Float64,rw,1);
m=200;
ngs = 3;          # No of gram-Schmidt orthogonalization
ifmod = false;    # If modified GS
println("nn=$nn; Matrix A($rw,$cl); m=$m")
if ifmod
  println("Using modified Gram-Schmidt orthogonolization")
else
  println("Using $ngs pass Gram-Schmidt orthogonolization")
end

Q,H = @timeit to "arnoldi" arnoldi(A,b,m,ngs,ifmod);

arnoldi_fac=norm(Q*H-A*Q[:,1:m]);
ortho=norm(Q'*Q-I);

println("Arnoldi Factorization condition satisfied to: $arnoldi_fac")
println("Q Matrix Orthogonality satisfied to: $ortho\n\n")


#evals = eigvals(H[1:m,1:m]);
#esort = sort(evals,by=abs,rev=true);

#eval_ref = eigvals(Matrix(A));
#esort2 = sort(eval_ref,rev=true);
#plot(esort2, linestyle="none",marker=".");
#
#plot(esort, linestyle="none",marker=".");
to



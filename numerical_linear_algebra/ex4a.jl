using Pkg
#Pkg.add("PyPlot")
#Pkg.add("LaTeXStrings")
#Pkg.add("PyCall")

using LinearAlgebra # for: eig, norm, etc
using MatrixDepot, Random,TimerOutputs,PyPlot,LaTeXStrings
using PyCall

include("arnoldi.jl")
include("GS_npass.jl")

println("Homework 1, Exercise 4a")

const to = TimerOutput();

# Create Matrix A
nn=100;
Random.seed!(0);        # To ensure I have consistent results across runs
A=randn(nn,nn); 
b=randn(nn);
rw,cl = size(A);
m=80;       # Final Krylov space dimension
ngs = 1;    # Single/double/Tripple/n Grahm Schmidt orthogonolizations
ifmod = false;
if ifmod
  println("Using modified Gram-Schmidt orthogonolization")
else
  println("Using $ngs pass Gram-Schmidt orthogonolization")
end

Q,H = @timeit to "arnoldi" arnoldi(A,b,m,ngs,ifmod);
should_be_zero1=norm(Q*H-A*Q[:,1:m]);
should_be_zero2=norm(Q'*Q-I);

println("Arnoldi Factorization condition satisfied to : $should_be_zero1")
println("Q Matrix Orthogonality satisfied to : $should_be_zero2\n\n")

lmbda = eigvals(H[1:m,1:m]);

PyCall.PyDict(matplotlib["rcParams"])["markers.fillstyle"] = "none"
#PyCall.PyDict(matplotlib["rcParams"])["lines.marker"] = "o"

pl = plot(real(lmbda),imag(lmbda), linestyle="none", marker="o");
ylabel(L"\lambda_{i}");
xlabel(L"\lambda_{r}");

to


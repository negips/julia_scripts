using Pkg
#Pkg.add("PyPlot")
#Pkg.add("LaTeXStrings")
#Pkg.add("PyCall")

using LinearAlgebra # for: eig, norm, etc
using MatrixDepot, Random,TimerOutputs,PyPlot,LaTeXStrings
using PyCall

include("arnoldi.jl")
include("GS_npass.jl")
include("exparnoldi.jl")

println("Homework 1, Exercise 4a")

close("all")

const to = TimerOutput();

# Create Matrix A
nn=100;
Random.seed!(0);        # To ensure I have consistent results across runs
U=randn(nn,nn);
diags = randn(nn);
b=rand(nn);
b[1] = -1;
b[2] = -1;
b[3] = -1;

B = diagm(0 => b);
B[1,2] = 1.;
B[2,3] = 1.;

A = inv(U)*B*U;

sort!(b);
e_exact = b;
ei = 0*b;

plot(b,ei,linestyle="none",marker=".",color="b")

rw,cl = size(A);
m=80;       # Final Krylov space dimension
ngs = 1;    # Single/double/Tripple/n Grahm Schmidt orthogonolizations
ifmod = false;
if ifmod
  println("Using modified Gram-Schmidt orthogonolization")
else
  println("Using $ngs pass Gram-Schmidt orthogonolization")
end

dt=0.5;
Q,H = @timeit to "arnoldi" exparnoldi(A,b,m,ngs,ifmod,dt);
should_be_zero1=norm(Q*H-A*Q[:,1:m]);
should_be_zero2=norm(Q'*Q-I);

println("Arnoldi Factorization condition satisfied to : $should_be_zero1")
println("Q Matrix Orthogonality satisfied to : $should_be_zero2\n\n")

ritz = eigvals(H[1:m,1:m]);
#sort!(lmbda)

PyCall.PyDict(matplotlib["rcParams"])["markers.fillstyle"] = "none"
#PyCall.PyDict(matplotlib["rcParams"])["lines.marker"] = "o"
#pl = plot(lmbda,linestyle="none", marker="o",color="red");

lmbda_r = log.(abs.(ritz))/dt;
lmbda_i = atan.(imag.(ritz),real.(ritz))/dt;

pl = plot(lmbda_r,lmbda_i, linestyle="none", marker="o",color="r");
#ylabel(L"\lambda_{i}");
#xlabel(L"\lambda_{r}");

to


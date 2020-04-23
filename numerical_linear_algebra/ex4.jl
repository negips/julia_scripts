using Pkg
#Pkg.add("PyPlot")
#Pkg.add("LaTeXStrings")
#Pkg.add("PyCall")

using LinearAlgebra # for: eig, norm, etc
using MatrixDepot, Random,TimerOutputs,PyPlot,LaTeXStrings
using PyCall

include("arnoldi.jl")
include("GS_npass.jl")
include("krylov_pm.jl")

println("Homework 1, Exercise 4a")

const to = TimerOutput();

# Create Matrix A
nn=10;
Random.seed!(0)
A=matrixdepot("wathen",nn,nn);
rw,cl = size(A);

b=randn(Float64,rw,1);
m=30;
ngs = 2;
ifmod = false;    # If modified GS
println("nn=$nn; Matrix A($rw,$cl); m=$m")
if ifmod
  println("Using modified Gram-Schmidt orthogonolization")
else
  println("Using $ngs pass Gram-Schmidt orthogonolization")
end

Km = @timeit to "krylov_pm" krylov_pm(A,b,m);

Q,H = @timeit to "arnoldi" arnoldi(A,b,m,ngs,ifmod);
should_be_zero1=norm(Q*H-A*Q[:,1:m]);
should_be_zero2=norm(Q'*Q-I);

println("Arnoldi Factorization condition satisfied to : $should_be_zero1")
println("Q Matrix Orthogonality satisfied to : $should_be_zero2\n\n")

lmbda_arnoldi = eigvals(H[1:m,1:m]);
lsort_a = sort(lmbda_arnoldi,by=abs,rev=true);

Am = Km'*A*Km;
Bm = Km'*Km;

lmbda_pmk = eigvals(Am,Bm);

lsort_p = sort(lmbda_pmk,by=abs,rev=true);

PyCall.PyDict(matplotlib["rcParams"])["markers.fillstyle"] = "none"
#PyCall.PyDict(matplotlib["rcParams"])["lines.marker"] = "o"

pla = plot(lsort_a, linestyle="-", marker="o",markersize=10,label="Arnoldi m=$m");
plp = plot(lsort_p, linestyle="--", marker="+",markersize=10,label="PM Kyrylov m=$m");
legend(fontsize=16)

ylabel(L"\lambda_{r}",fontsize=20);
xlabel("""Krylov space "m" """,fontsize=16);

#savefig("pm_arnoldi.eps")
to


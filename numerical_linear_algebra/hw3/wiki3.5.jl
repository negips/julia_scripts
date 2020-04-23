#!/usr/bin/julia
# Source code for Homework 3

println("Homework 3, Wiki question 3.5")

#import Pkg
#Pkg.add("PyPlot")
#Pkg.add("PyCall")
#Pkg.add("Blink")

using Blink,TimerOutputs
using PyPlot,PyCall,LinearAlgebra,Random

include("ShiftedQR_err.jl")
include("QRsimple.jl")
include("QRsimple_err.jl")
include("QRfac.jl")
include("GS_modified.jl")
include("GS_npass.jl")
include("geterr.jl")

close("all")      # close all open figures

const to = TimerOutput();

global A
m=20;
Random.seed!(0)

A = rand(Float64,m,m);
A = A + A';

B = copy(A);
eigs = eigvals(A);
eigs = sort(eigs);

pl = plot(eigs, marker="o", linestyle="none", label="\$\\lambda\$");
ylabel(L"\lambda",fontsize="14")
# xlabel(L"\mathcal{R}",fontsize="14")
legend(fontsize=14)


ifmod = false;    # if we use modified GS
ngs = 2;          # No of GS passes
errtol = 1.0e-14; # Tolerance for off-diagonal term
nstep = 10;       # No of iteration steps. 
                  # Keep negative if we want to iterate to specified tolerance.
neig=m;   

# Perform simple QR
A,errhist,eigerr1 = @timeit to "QR" QRsimple_err(A,errtol,ifmod,ngs,eigs);
evals1 = sort(diag(A));

plot(evals1, marker="s", linestyle="none", label="QR simple")

ifrayleigh=true
shift = B[m,m];
H,evals,errhist2,eigerr2 = @timeit to "RayleighQR" ShiftedQR_err(B,neig,shift,ifrayleigh,nstep,errtol,eigs);
#semilogy(eigerr2)

evals2 = sort(diag(H));

plot(evals2, marker=".", linestyle="none", label="Rayleigh QR")
legend()

figure()
semilogy(errhist2)
#semilogy(eigerr2)

#fname="ex3d_h21.png";
#savefig(fname)

to

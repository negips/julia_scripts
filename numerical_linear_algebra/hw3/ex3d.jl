#!/usr/bin/julia
# Source code for Homework 3

println("Homework 3, Exercise 3d")

#import Pkg
#Pkg.add("PyPlot")
#Pkg.add("PyCall")
#Pkg.add("Blink")

using Blink,TimerOutputs
using PyPlot,PyCall,LinearAlgebra

include("ShiftedQR.jl")
include("QRsimple.jl")
include("QRfac.jl")
include("GS_modified.jl")
include("GS_npass.jl")
include("geterr.jl")

close("all")      # close all open figures

const to = TimerOutput();

ifmod = false;    # if we use modified GS
ngs = 1;          # No of GS passes

errtol = 1.0e-10; # error tolerance for sub diagonal elements

alpha = 1;
e_array = [0.4 0.1 0.01 1.0e-3 1.0e-4 1.0e-5 1.0e-6 1.0e-7 1.0e-8 1.0e-9 1.0e-10 0];
ne = length(e_array);
e_array = reshape(e_array,ne,1);
nstep = 1;       # No of iteration steps. 
                  # Keep negative if we want to iterate to specified tolerance.
errtol = 1.0e-12; # Tolerance for off-diagonal term
neig=1;   

shifts = [0. 1.];
nshifts = length(shifts);

global H
H = [3. 2.; 0 1];

evals_array = zeros(ne,nshifts);
h21 = zeros(ne,nshifts);

for ii in 1:nshifts
  global shift
  shift = shifts[ii];

  for jj in 1:ne
    global U
  
    e = e_array[jj];
    H[2,1] = e;
  
    U,evals = @timeit to "ShiftedQR" ShiftedQR(H,neig,shift,nstep,errtol);
    evals_array[jj,ii] = evals[1]; 
    h21[jj,ii] = U[2,1];
  end

  pl = loglog(e_array,h21[:,ii], marker="o", label="Shift=$shift");
  ylabel(L"\bar{h}_{2,1}",fontsize="14")
  xlabel(L"\epsilon",fontsize="14")
  legend(fontsize=14)
 
end
fname="ex3d_h21.png";
savefig(fname)

to

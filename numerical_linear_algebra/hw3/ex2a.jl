#!/usr/bin/julia
# Source code for Homework 3

println("Homework 3, Exercise 2a")

#import Pkg
#Pkg.add("PyPlot")
#Pkg.add("PyCall")
#Pkg.add("Blink")

using Blink
using PyPlot,PyCall,LinearAlgebra

include("QRfac.jl")
include("GS_npass.jl")
include("GS_modified.jl")
include("geterr.jl")
include("alpha_example.jl")
include("QRsimple.jl")

close("all")      # close all open figures

ifmod = false;    # if we use modified GS
ngs = 1;          # No of GS passes

errtol = 1.0e-10; # error tolerance for sub diagonal elements

powers = [i for i in -1:6];
alpha = [0.5 1 2 5 10 25 50 100 250 500 1000 2500 5000 10000 25000 50000 100000 250000 500000];
#alpha = 5*(10.0.^powers);
#alpha = powers;
nalpha = length(alpha);
alpha = reshape(alpha,nalpha,1);
itersalpha = zeros(nalpha,1);
iters_Theory = zeros(nalpha,1);
Gammafac = zeros(nalpha,1);
Cfac = zeros(nalpha,1);

Niters = 10000;

for ia in 1:nalpha

  alphai = alpha[ia];

# Create Matrix A
  global A
  A=alpha_example(alphai,20);
  evals = eigvals(A);
  se=sort(abs.(evals));
  global se
  gamma = 0.;
  global gamma
  neig=length(se);
  for i in 1:neig-1
    rmax2 = maximum(se[i]./se[i+1:neig]);
    gamma = max(gamma,rmax2);
  end
  iters_Theory[ia] = log(errtol)/log(gamma);
  Gammafac[ia]=gamma;

  A,errhist = QRsimple(A,errtol,ifmod,ngs);
  nitrs = length(errhist)
  itersalpha[ia]=nitrs;

  Cfac[ia] = exp((iters_Theory[ia]-itersalpha[ia])*log(gamma))
  
#  layout = ["yaxis" => ["type"=>"log", "autorange" => true] ];
#  # pl = plot(err_hist, ["layout" => layout])
#  pl = semilogy(errhist, label="$alphai");
#  ylabel(L"|\epsilon|",fontsize="14")
#  xlabel("Iterations",fontsize="14")
#  legend()

end

h2 =figure();
pl3 = semilogx(alpha,iters_Theory, linestyle="--", marker="o", label="Prediction");
pl2 = semilogx(alpha,itersalpha, linestyle="-", marker="o", color="r", label="Computed");

xlabel(L"\alpha", fontsize="14");
ylabel("Iterations",fontsize="14");
fname="ex2a_iterations.png"
legend()

savefig(fname)

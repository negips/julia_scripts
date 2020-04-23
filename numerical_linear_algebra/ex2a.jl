#!/usr/bin/julia
# Source code for Homework 1

println("Homework 1, Exercise 2a")

#import Pkg
#Pkg.add("PyPlot")
#Pkg.add("PyCall")
#Pkg.add("Blink")

using Blink
using PyPlot,PyCall,LinearAlgebra


# Create Matrix A
A = [1 2 3; 2 2 2; 3 2 9];

evals = eigvals(A);

mu_ref = maximum(abs.(evals));

# Initial vector
x0 = [1,1,1];

# Power Method

v = copy(x0);
vnorm = sqrt(v'*v);
v = v./vnorm;
vp = copy(v); # retain copy of last vector
err_norm = 1.;     # Initialization

Niters = 40;       # No of iterations
err_hist = Vector{Float64}(undef,Niters);     # Just initialize

for i in 1:Niters
   global v
   global vp
   global err_hist

#   eps = 1.0e-20;

   v = A*v
   vnorm = sqrt(v'*v);
   v = v./vnorm;
#   err_norm = sqrt((vp - v)'*(vp - v));
   mu = v'*A*v./(v'*v);
   err_norm = max((mu_ref - mu),eps());

   err_hist[i] = err_norm;
   vp = copy(v);
end

#deleteat!(err_hist,1)

j = [i for i in 1:Niters];

layout = ["yaxis" => ["type"=>"log", "autorange" => true] ];
# pl = plot(err_hist, ["layout" => layout])
pl = semilogy(err_hist)
ylabel(L"|\epsilon|")
xlabel("Iterations")

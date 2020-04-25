#!/usr/bin/julia
# Source code for Exercise 2 of Homework 1

println("Matrix Norm of a Jordan form Matrix")

import Pkg
#Pkg.add("PyPlot")
#Pkg.add("Blink")

using Blink
using PyPlot
using LinearAlgebra


# Create Matrix A
v = [-1.0, -1.0, -0.1, -0.1, -0.01, -0.01];

n = length(v);
A = zeros(Float64,n,n);

for i in 1:n
  A[i,i] = v[i];
end

#epsilon = 1.0e-1;
#A[1,2] = epsilon;

epsilon = 1.00e-0;
A[3,4] = epsilon;

epsilon = 0.1e-0;
A[5,6] = epsilon;


# v = copy(x0);
vnorm = v'*v;
v = v./vnorm;
vp = copy(v); # retain copy of last vector
err_norm = 1.;     # Initialization

t = range(1.0e-6, step=1.0e-2, stop=500.0);

Niters = length(t);

err_hist = Vector{Float64}(undef,Niters);     # Just initialize

for i in 1:Niters
   global v
   global vp
   global err_hist

   v = exp(A*t[i]);
   vnorm = norm(v);
   err_hist[i] = vnorm;
end

#deleteat!(err_hist,1)

# Not sure what this is called. Loading a module?
# using Plotly

j = [i for i in 1:Niters];

layout = ["yaxis" => ["type"=>"log", "autorange" => true] ];
# pl = plot(err_hist, ["layout" => layout])
pl = plot(t,err_hist)
ylabel("|exp(At)|");
xlabel("time");

println("Done.")


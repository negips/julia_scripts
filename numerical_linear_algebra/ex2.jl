#!/usr/bin/julia
# Source code for Exercise 2 of Homework 1

println("Homework 1, Exercise 2")

import Pkg
Pkg.add("PyPlot")
Pkg.add("Blink")

using Blink
using PyPlot


# Create Matrix A
A = [1 2 3; 2 2 2; 3 2 9];

# Initial vector
x0 = [1,1,1]

# Power Method

v = copy(x0);
vnorm = v'*v;
v = v./vnorm;
vp = copy(v); # retain copy of last vector
err_norm = 1.;     # Initialization

Niters = 40;       # No of iterations
err_hist = Vector{Float64}(undef,Niters);     # Just initialize

for i in 1:Niters
   global v
   global vp
   global err_hist

   v = A*v
   vnorm = sqrt(v'*v);
   v = v./vnorm;
   err_norm = sqrt((vp - v)'*(vp - v));
   err_hist[i] = err_norm;
   vp = copy(v);
end

#deleteat!(err_hist,1)

# Not sure what this is called. Loading a module?
# using Plotly

j = [i for i in 1:Niters];

layout = ["yaxis" => ["type"=>"log", "autorange" => true] ];
# pl = plot(err_hist, ["layout" => layout])
pl = semilogy(err_hist)
ylabel("|err|")
xlabel("Iterations")

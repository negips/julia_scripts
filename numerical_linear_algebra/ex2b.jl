#!/usr/bin/julia
# Source code for Homework 1

println("Homework 1, Exercise 2b")

import Pkg
#Pkg.add("PyPlot")
#Pkg.add("Blink")
#Pkg.add("LinearAlgebra")
#Pkg.add("IJulia")

using IJulia
#IJulia.installkernel("Julia nodeps", "--depwarn=no")

using Blink
using PyPlot
using LinearAlgebra

# Create Matrix A
A = [1 2 3; 2 2 2; 3 2 9];

# Find Eigevnalues for reference     
evals = eigvals(A);

# Initial vector
x0 = [1,1,1]

# Rayleigh Quotient Iteration

v = copy(x0);
vnorm = sqrt(v'*v);
v = v./vnorm;
vp = copy(v); # retain copy of last vector

Niters = 10;       # No of iterations
err_hist = Vector{Float64}(undef,Niters);     # Just initialize
mu_hist  = Vector{Float64}(undef,Niters);     # Just initialize

for i in 1:Niters
   global v
   global vp
   global err_hist
   global mu_hist
   global mu
   global vn
   local epsilon

   epsilon=1.0e-20;     # Just something below machine precision

   eye = Matrix{Float64}(I,3,3);
   mu = (v'*A*v)/(v'*v);

   mu_hist[i] = mu;

   B = inv(A - mu*eye);
   w = B*v;
   wnorm = sqrt(w'*w);
   v = w/wnorm;

#  Eigenvector change  
   err_norm1 = sqrt((vp-v)'*(vp-v));
   err_norm2 = sqrt((vp+v)'*(vp+v));

   if err_norm1<err_norm2
     err_norm = err_norm1+epsilon
   else
     err_norm = err_norm2+epsilon
   end

   err_hist[i] = err_norm;

   vp = copy(v);

end

j = [i for i in 1:Niters];

println("Eigenvalue estimate after $Niters iterations: $mu")
println("""Eigenvalues Calculated using "eigvals" function:= $evals """)

index = argmin(abs.(evals .- mu));
mu_error = abs.(mu_hist .- evals[index]) .+ 1.0e-20;

pl = semilogy(mu_error);
ylabel("|err|");
xlabel("Iterations");



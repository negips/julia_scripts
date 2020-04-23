"""

Example code to verify that the eigenvectors are stationary points of the 
Rayleigh quotient in the symmetric case

"""

# Load the necessary modules
using LinearAlgebra
using Plots
using MatrixDepot, Random, Arpack
using MAT

N = 100;

# Restart the seed
Random.seed!(0);

A = rand(N,N);
A = (A+A'); # Forming a symmetric matrix

eigval,eigvec = eigen(A);

for i in 1:N
    x = eigvec[:,i];
    gradx = (2*A*x*(x'*x)-2*(x'*A*x)*x)/((x'*x).^2);
    println("norm(gradx)= ",norm(gradx))
end

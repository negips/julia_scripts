# Load the necessary modules
using LinearAlgebra
using Plots
using MatrixDepot, Random, Arpack
using MAT


########################################
#           EXERCISE 3                 #
########################################
include("arnoldi.jl")

# Restart the seed
Random.seed!(0);

# Set up the matrices
nn=10;
A = matrixdepot("wathen",nn,nn);
b = randn(3*nn*nn+2*nn+2*nn+1);

# Define the number of iterations
m = 100;
Q,H= @time arnoldi(A,b,m);

# Plot the eigenvalues for A and H
λ, ϕ  = eigs(A, nev = 120);
p2 = scatter(λ,marker=2)
       
λH, ϕH  = eigen(H[1:m,1:m]);
sort!(λH, by=abs, rev=true)
p2 = scatter!(p2,real(λH),marker=3);
display(plot(p2))

should_be_zero1=norm(Q*H-A*Q[:,1:m])
orth=norm(Q'*Q-I)

########################################
#           END EXERCISE 3             #
########################################

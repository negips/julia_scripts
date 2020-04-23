########################################
#           EXERCISE 4                 #
########################################

# Load the necessary modules
using LinearAlgebra
using Plots
using MatrixDepot, Random, Arpack
using MAT


# Restart the seed
Random.seed!(0);

# Set up the matrices
nn=10;
A = matrixdepot("wathen",nn,nn);
b = randn(3*nn*nn+2*nn+2*nn+1);

# Do in a loop to generate the figure
maxm  = 10;
λmaxK = zeros(maxm,1);
λmaxA = zeros(maxm,1);

for m in 1:maxm
    println("m ",m)
    # Generate the Km matrix
    Km = zeros(length(b),m);
    
    for i in 1:m
        Km[:,i] = (A^(i-1)*b)/norm(A^(i-1)*b,2); 
    end
    
    # Solve the generalized eigenvalue problem Ax = λBx
    λ4,ϕ4 = eigen(Km'*A*Km,Km'*Km);
    sort!(real(λ4), by=abs, rev=true)
    λmaxK[m] = maximum(real(λ4));

    # Done the prior method
#    println("Done the prior method, will do Arnoldi")
#    readline(stdin)
    
    
    # Do the Arnoldi for comparison
    Q,H= arnoldi(A,b,m);
    λH, ϕH  = eigen(H[1:m,1:m]);
    sort!(λH,by=abs, rev=true)
    λmaxA[m] = maximum(abs.(λH));
    
end

p3 = plot(collect(1:1:maxm),λmaxA,marker=3,label="Arnoldi")
p3 = scatter!(p3,collect(1:1:maxm),λmaxK,marker=2,label="Approx in (2)")
display(plot(p3))

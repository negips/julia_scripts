# Code for Julia program

println("Julia script Hw1")

# Load the necessary modules
using LinearAlgebra
using Plots
using MatrixDepot, Random, Arpack
using MAT


########################################
#           EXERCISE 2                 #
########################################

include("functionsEx2.jl")

# Matrix from HW1 exercise 2
A    = [1 2 4;2 2 2; 3 2 9];

x0   = [1;1;1];     # Initial guess
maxk = 20;          # Max iterations

Avals, Avecs = eigen(A);
maxAvals = maximum(Avals);
sort!(Avals);

maxEigPM,vPM,itEigPM=powerIt(A,x0,maxk);
itErrPM = norm.(itEigPM.-maxAvals,2);

# There is a problem with the error being exactly 0 at some points
# presenting an issue with the log plot. The next lines are done
# to fix the error to numerical epsilon if the eigErr is too low

for j=1:length(itErrPM)
    itErrPM[j]=max(itErrPM[j],eps())
end

maxEigRQI,vRQI,itEigRQI=RQit(A,x0,maxk);
itErrRQI = abs.(itEigRQI.-maxAvals);

# There is a problem with the error being exactly 0 at some points
# presenting an issue with the log plot. The next lines are done
# to fix the error to numerical epsilon if the eigErr is too low

for j=1:length(itErrRQI)
    itErrRQI[j]=max(itErrRQI[j],eps())
end

gr()
p1=plot(collect(1:1:maxk),itErrPM,label="eigErr in PM",yaxis=:log10,marker=2)
p1=plot!(p1,collect(1:1:maxk),itErrRQI,label="eigErr in RQI",yaxis=:log10,marker=2)

xlabel!("Iterations")
ylabel!("\$ \\lambda-\\lambda^k \$")
ylims!((1e-16,1e1))

########################################
#           END EXERCISE 2             #
########################################


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
sort!(λH, rev=true)
p2 = scatter!(p2,λH,marker=2);
#p2 = scatter(λH,marker=2);
#display(plot(p1,p2,layout=(2,1)))

should_be_zero1=norm(Q*H-A*Q[:,1:m])
orth=norm(Q'*Q-I)

########################################
#           END EXERCISE 3             #
########################################


########################################
#           EXERCISE 4                 #
########################################

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
    sort!(real(λ4), rev=true)
    λmaxK[m] = maximum(real(λ4));

    # Done the prior method
#    println("Done the prior method, will do Arnoldi")
#    readline(stdin)
    
    
    # Do the Arnoldi for comparison
    Q,H= arnoldi(A,b,m);
    λH, ϕH  = eigen(H[1:m,1:m]);
    sort!(λH, rev=true)
    λmaxA[m] = maximum(λH);
    
end

p3 = plot(collect(1:1:maxm),λmaxA,marker=3,label="Arnoldi")
p3 = plot!(p3,collect(1:1:maxm),λmaxK,marker=2,label="Approx in (2)")

#display(plot(p1,p2,p3,layout=(3,1)))

########################################
#           END EXERCISE 4             #
########################################


########################################
#           EXERCISE 6                 #
########################################

include("exercise6.jl")

println("Julia script Hw1 DONE!")

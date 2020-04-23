"""

Example for the classical vs. the modified GS

"""

# Load the necessary modules
using LinearAlgebra
using Plots
using MatrixDepot, Random, Arpack
using MAT

function CGS(A,w)

    # Algorithm for the classical Gram-Schmidt

    # Input
    # A: Matrix on which to apply the GS algorithm
    # w: vector to orthogonalise with respect to the basis
    
    # Output
    # A_GS: Is the matrix on which GS has been applied
    # the output is not normalised
    
    n,m  = size(A);
    A_GS = zeros(n,m+1);

#    orth_err = norm(A'A-Matrix{Float64}(I, m, m));
        
#    if orth_err >= 1e5
#        println("The given basis is not orthongonal A'A neq. I")
#        return
#    end

    y = w;
    h = zeros(m,1); # Initialize h outside the loop
    for i in 1:m
        h[i] = A[:,i]'*w;
        y = y-A[:,i].*h[i];
    end
    β = norm(y,2);
    y = y/β;

    A_GS = hcat(A,y);
    return A_GS
    
end

function MGS(A,w)

    # Algorithm for the modified Gram-Schmidt

    # Input
    # A: Matrix on which to apply the GS algorithm
    # w: vector to orthogonalise with respect to the basis
    
    # Output
    # A_GS: Is the matrix on which GS has been applied
    # the output is not normalised
    
    n,m  = size(A);
    A_GS = zeros(n,m+1);

    orth_err = norm(A'A-I);
        
    if orth_err >= 1e5
        println("The given basis is not orthongonal A'A neq. I")
        return
    end

    y = w;
    h = zeros(m,1); # Initialize h outside the loop
    for i in 1:m
        h[i,1] = A[:,i]'*y;
        y = y-A[:,i]*h[i];
    end
    β = norm(y,2); 
    y = y/β;

    A_GS = hcat(A,y);
    return A_GS
    
end

epsi = 1e-12; # Noise to the initial matrix
Q = [1.0 0;epsi 1.0;0 0;0 0]; # Initial orthogonal matrix

#Q = Q.+epsi*rand(4,1);

w = [1.0; 0; 0; epsi]; # Vector to be added to the base
#w = w.+epsi*rand(Float64,1,1);

Q_CGS = CGS(Q,w);
Q_MGS = MGS(Q,w);

should_be_zeroCGS = norm(Q_CGS'Q_CGS-I);
println("Q,w,Q_CGS",Q,w,Q_CGS)
println("CGS is 0? -> ",should_be_zeroCGS)

should_be_zeroMGS = norm(Q_MGS'Q_MGS-I);
println("Q,w,Q_MGS",Q,w,Q_MGS)
println("CGS is 0? -> ",should_be_zeroMGS)






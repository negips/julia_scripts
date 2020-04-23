"""
Does Arnoldi with a shift

Q,H=arnoldi(A,b,m)
A simple implementation of the Arnoldi method.
The algorithm will return an Arnoldi "factorization":
Q*H[1:m+1,1:m]-A*Q[:,1:m]=0
where Q is an orthogonal basis of the Krylov subspace
and H a Hessenberg matrix.

Example:
```julia-repl
A=randn(100,100); b=randn(100,1);
m=10;
Q,H=arnoldi(A,b,m);
println("should_be_zero1=",norm(Q*H-A*Q[:,1:m]))
println("should_be_zero2=",norm(Q'*Q-I))
```

"""
function arnoldi_shift(A,b,mu,m)

    n=length(b);
    Q=zeros(Complex{Float64},n,m+1);
    H=zeros(Complex{Float64},m+1,m);
    Q[:,1]=b/norm(b);

    t0=0;
    t1=0;

    for k=1:m
        
        w=(A-mu*I)\Q[:,k]; # Matrix-vector product with last element
        # Orthogonalize w against columns of replace this with a orthogonalization
        Q[:,1:k+1],β,h=DGS(Q[:,1:k],w);
        #Put Gram-Schmidt coefficients into H
        H[1:(k+1),k]=[h;β];
        # normalize
        # OBS the normalization is done within GS variant
        # Plot the eigenvalues for A and H


    end
    return Q,H
end

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

    orth_err = norm(A'A-Matrix{Float64}(I, m, m));
        
    if orth_err >= 1e5
        println("The given basis is not orthongonal A'A neq. I")
        return
    end

    y = w;
    h = zeros(m,1); # Initialize h outside the loop
    for i in 1:m
        h[i] = A[:,i]'*w;
        y = y-A[:,i]*h[i];
    end
    β = norm(y,2); 
    y = y/β;

    A_GS = hcat(A,y);
    return A_GS,β,h
    
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
        h[i] = A[:,i]'*y;
        y = y-A[:,i]*h[i];
    end
    β = norm(y,2); 
    y = y/β;

    A_GS = hcat(A,y);
    return A_GS,β,h
    
end

function DGS(A,w)

    # Algorithm for the double Gram-Schmidt

    # Input
    # A: Matrix on which to apply the GS algorithm
    # w: vector to orthogonalise with respect to the basis
    
    # Output
    # A_GS: Is the matrix on which GS has been applied
    # the output is not normalised
   
    n,m  = size(A);
    A_GS = zeros(Complex{Float64},n,m+1);

    orth_err = norm(A'A-I);
        
    if orth_err >= 1e5
        println("The given basis is not orthongonal A'A neq. I")
        return
    end

    y = w;
    h = zeros(Complex{Float64},m); # Initialize h outside the loop
    g = zeros(Complex{Float64},m); # Initialize h outside the loop

    for i in 1:m
        h[i] = A[:,i]'*w;
        y = y-A[:,i]*h[i];
    end
    w = y;
    for i in 1:m
        g[i] = A[:,i]'*w;
        y = y-A[:,i]*g[i];
        h[i] = h[i]+g[i];
    end
    
    β = norm(y,2); 
    y = y/β;
   
    A_GS = hcat(A,y);
    
    return A_GS,β,h
    
end

function TGS(A,w)

    # Algorithm for the triple Gram-Schmidt

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
    g = zeros(m,1);
    f = zeros(m,1);
    for i in 1:m
        h[i] = A[:,i]'*w;
        y = y-A[:,i]*h[i];
        g[i] = A[:,i]'*y;
        y = y-A[:,i]*g[i];
        f[i] = A[:,i]'*y;
        y = y-A[:,i]*f[i];
        h[i] = h[i]+g[i]+f[i];
    end
    β = norm(y,2); 
    y = y/β;

    A_GS = hcat(A,y);
    return A_GS,β,h
    
end

using LinearAlgebra

function naive_hessenberg_red(A)
    # A naive inefficient way of carrying out
    # the hessenberg reduction
    #
    n=size(A,1)
    for k=1:n-2
        x=A[k+1:n,k];
        xnorm = norm(x);    
        y=zeros(n-k,1);
        y[1] = xnorm;
        alpha = 1;    
        z=x - alpha*y; 
        u=z/norm(z);
        P1=I-2*u*u'
        P0=Matrix{Float64}(I, k, k)
        P=[P0 zeros(k,n-k); zeros(n-k,k) P1]  # Equation (2.5) in lecture notes
        PA=P*A
        # PA should have zeros below the second off-diagonal in column k as described on page 6 in lecture notes.
        A=PA*P'
        # A should have the same zero-structure as PA
    end
    return H=A # should be a hessenberg matrix with same eigenvalues as input A
end

#A=rand(5,5)
#H=naive_hessenberg_red(A);

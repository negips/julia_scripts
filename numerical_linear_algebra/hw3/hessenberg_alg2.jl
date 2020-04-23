using LinearAlgebra

function hessenberg_alg2(A)
    # Algorithm 2 for hessenberg reduction
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
#       Compute [Pk]*[A]
        A[k+1:n,k:n]=A[k+1:n,k:n]- 2*u*(u'*A[k+1:n,k:n]);
#       Compute [Pk]*[A]*[Pk]^T
        A[1:n,k+1:n]=A[1:n,k+1:n]- 2*(A[1:n,k+1:n]*u)*u';
    end
    return A # should be a hessenberg matrix with same eigenvalues as input A
end

#A=rand(5,5)
#H=naive_hessenberg_red(A);

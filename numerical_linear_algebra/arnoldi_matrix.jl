"""
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
function arnoldi_matrix(A,m,ngs,ifmod)

    b = A[:,1];  
    n=length(b);
    Q=zeros(n,m);
    H=zeros(m,m);
    Q[:,1]=b/norm(b);

    t0=0;
    t1=0;

    if ifmod
      Q,H=GS_modified(A,m);
      return Q,H
    end  

    for k=2:m
       w=A[:,k];  # Orthogonalize columns of A
#      Orthogonalize w against columns of Q
       h,β,worth=GS_npass(Q,w,k-1,ngs);
#      Put Gram-Schmidt coefficients into H
       H[1:k,k]=[h;β];
#      normalize
       Q[:,k]=worth/β;
    end

    return Q,H
end

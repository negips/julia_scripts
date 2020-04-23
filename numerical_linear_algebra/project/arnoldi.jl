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
function arnoldi(A,b,m,ngs,ifmod)

    n=length(b);
    if m>n
      m=n
    end  

    Q=zeros(n,m+1);
    H=zeros(m+1,m);
    Q[:,1]=b/norm(b);

    t0=0;
    t1=0;

    for k=1:m
        w=A*Q[:,k]; # Matrix-vector product with last element
#       Orthogonalize w against columns of Q
        if ifmod
          h,β,worth=GS_modified(Q,w,k);
        else
          h,β,worth=GS_npass(Q,w,k,ngs);
        end
#       Put Gram-Schmidt coefficients into H
        H[1:(k+1),k]=[h;β];
#       normalize
        Q[:,k+1]=worth/β;
    end

    return Q,H
end

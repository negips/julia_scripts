function mygmres(A,b,ngs,ifmod,tol,maxit,xexact,ifcomplex)

    m = maxit;      # Assuming I don't take more than 100 iterations

    n=length(b);
    if m>n
      m=n
    end

    if (ifcomplex)
      Q=zeros(Complex,n,m+1);
      H=zeros(Complex,m+1,m);
    else
      Q=zeros(Float64,n,m+1);
      H=zeros(Float64,m+1,m);
    end

    residuals = zeros(Float64,m+1,1);
    xerr      = zeros(Float64,m+1,1);

    bnorm = norm(b); 
    Q[:,1]=b/bnorm;
    xnew = 0*b;  

    tol=abs(tol);
    k = 0;
    resid = 1.
    while k<m && resid>tol
#        global z
        k = k+1;
        w=A*Q[:,k]; # Matrix-vector product with last element

#       Orthogonalize w against columns of Q
        if ifmod
          h,β,worth=GS_modified(Q,w,k,ifcomplex);
        else
          h,β,worth=GS_npass(Q,w,k,ngs,ifcomplex);
        end
#       Put Gram-Schmidt coefficients into H
        H[1:(k+1),k]=[h;β];
#       normalize
        Q[:,k+1]=worth/β;

#       Find least-squares approximation to ||Hk*z - |b|e1||
        be1 = zeros(k+1,1);
        be1[1] = bnorm;
        Qh,Rh = qr(H[1:k+1,1:k]);
        Qh = Qh[1:k+1,1:k];
        bq = Qh'*be1;
        z = inv(Rh)*bq;
#        resid = norm(Qh*Qh'*be1 - be1);
        resid = norm(H[1:k+1,1:k]*z - be1); 
#        println("Residual norm after iteration $k = $resid")
        residuals[k]=resid;
        xnew = Q[:,1:k]*z;
        xerr[k] = norm(xexact-xnew);

    end

    return xnew,Q[:,1:k+1],H[1:k+1,1:k],residuals[1:k],xerr[1:k]
end

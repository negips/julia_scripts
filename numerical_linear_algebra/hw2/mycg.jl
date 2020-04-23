function mycg(A,b,tol,maxit,xexact,ifcgn,ifcomplex)

    m = maxit;      # Assuming I don't take more than 100 iterations

    n=length(b);
    if m>n
      m=n
    end 

    Ac = A';
    bold = b;

    if (ifcgn)          # If CG normal
      b = Ac*b;
    end  

    if ifcomplex
      X           = zeros(Complex,n,m);       # CG iterates
      P           = zeros(Complex,n,m+1);     # Search directions
      T           = zeros(Complex,m+1,m);     # Tri-diagonal matrix
      R           = zeros(Complex,n,m+1);     # Residual rectors
    else
      X           = zeros(Float64,n,m);       # CG iterates
      P           = zeros(Float64,n,m+1);     # Search directions
      T           = zeros(Float64,m+1,m);     # Tri-diagonal matrix
      R           = zeros(Float64,n,m+1);     # Residual rectors
    end  

    residuals   = zeros(Float64,m,1); 
    xerr        = zeros(Float64,m,1);

    
    xk     = 0*b;         # Initial vecotor is zero
    rk     = b;             # starting with x0 = 0;
    R[:,1] = rk;          
    pk     = rk;            # Initial search direction
    P[:,1] = rk;
    Ap     = xk;              # Just to make it global


    tol=abs(tol);
    k = 0;
    resid = 1.
    rr = rk'*rk;
    while k<m && resid>tol
        k         = k+1;
       
        if (ifcgn)
          Ap      = A*pk;
          Ap      = Ac*Ap;
        else
          Ap      = A*pk;
        end
        αk        = rr/(pk'*Ap);
        xk        = xk + αk.*pk;
        X[:,k]    = xk;

        rk        = rk - αk.*Ap;
        R[:,k+1]  = rk; 

        rnew      = rk'*rk;
        
        βk        = rnew/rr;
        pk        = rk + βk.*pk; 
        P[:,k+1]  = pk;

        rr              = rnew;
        resid           = norm(rk);
        residuals[k]    = resid;
#        if (ifcgn)
#          resid         = norm(A*xk-bold);
#          residuals[k]  = resid;
#        end  

        xerr[k]         = norm(xexact-xk);
    end          

    return xk,X[:,1:k],P[:,1:k+1],R[:,1:k+1],residuals[1:k],xerr[1:k]
end




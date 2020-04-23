using LinearAlgebra
function QRsimple_err(A,errtol,ifmod,ngs,eigs)
# Input:
#  A - input Matrix

  Niters     = 50^4;
  errhist    = zeros(Niters,1);
  errlast    = geterr(A);
  errhist[1] = errlast;
  eigerr     = zeros(Niters,1);

# I assume the eigenvalues are approximated by the diagonal elements.
# Even though the iterate matrix is not upper triangular

  eigs = sort(eigs);

  i=0;
  while errlast > errtol
    i=i+1;
#   Perform QR factorization
    Q,R = QRfac(A,ifmod,ngs);
    A = R*Q;
    approx_eig = copy(sort(diag(A)));
    eigerr[i] = maximum(abs.(approx_eig - eigs));  

    errlast = geterr(A);  
    errhist[i]=errlast;
  end

  return A,errhist[1:i],eigerr[1:i]

end


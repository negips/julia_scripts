using LinearAlgebra
function QRsimple(A,errtol,ifmod,ngs)
# Input:
#  A - input Matrix

  Niters     = 10^4;
  errhist    = zeros(Niters+1,1);
  errlast    = geterr(A);
  errhist[1] = errlast;

  i=1;
  while errlast > errtol
    i=i+1;
#    global A
#   Perform QR factorization
    Q,R = QRfac(A,ifmod,ngs);
    A = R*Q;
    errlast = geterr(A);  
    errhist[i]=errlast;
  end

  return A,errhist[1:i]

end


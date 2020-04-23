using LinearAlgebra
function QRfac(A,ifmod,ngs)
# Input:
#  A - input Matrix

  (n,m) = size(A);
  Q=zeros(n,m);
  R=zeros(m,m);
  b = A[:,1];
  Q[:,1]=b/norm(b);

  for k=1:m
      w=A[:,k]; # Matrix-vector product with last element
#     Orthogonalize w against columns of Q
      if ifmod
        h,β,worth=GS_modified(Q,w,k-1);
      else
        h,β,worth=GS_npass(Q,w,k-1,ngs);
      end
#     Put Gram-Schmidt coefficients into R
      R[1:k,k]=[h;β];
#     normalize
      Q[:,k]=worth/β;
  end

  return Q,R

end


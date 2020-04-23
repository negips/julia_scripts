using LinearAlgebra
function ShiftedQR(A,neig,shift,nstep,errtol)
# Input:
#  A - input Hessenberg Matrix
#  neig - No of eigenvalues to converge
#  (We start from the last row of H)
#  shift - shift applied to the Hessenberg Matrix
#  nstep - How many iterations to perform.
#  (Use negative value if we wish to converge to tolerance)
#  errtol - specified tolerance of eigenvalues.

  H = copy(A);
  nn = size(H);
  n = nn[1];
  global evals
  evals = zeros(neig,1);

  rend = n-neig+1;
# Going in Reverse order

  global i
  i = 0;
  for r in n:-1:rend
    global I
    global H2
    global k
    global offdiag  

    k=0;
    msize = r;          # Matrix size
    H2 = H[1:r,1:r];
    eye = Matrix{Float64}(I,2,2);     # Identity matrix
    offdiag = H2[r,r-1]; 

    while offdiag > errtol
      k=k+1;

      H2 = H2 - eye*shift;
      ifmod = false;          # if modified GS
      ngs = 1;                # no of GS passes.

      Q,R = QRfac(H2,ifmod,ngs);
      H2 = R*Q + eye*shift;
      offdiag = abs(H2[r,r-1]);
      
      if nstep>0 && k==nstep
        offdiag=0.
      end

    end

    i=i+1;
    evals[i]=H2[r,r];  
    H[1:r,1:r] = H2;  
  end

  return H,evals

end


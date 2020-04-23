using LinearAlgebra
function ShiftedQR_err(A,neig,shift,ifrayleigh,nstep,errtol,eigs)
# Input:
#  A - input Hessenberg Matrix
#  neig - No of eigenvalues to converge
#  (We start from the last row of H)
#  shift - shift applied to the Hessenberg Matrix
#  nstep - How many iterations to perform.
#  (Use negative value if we wish to converge to tolerance)
#  errtol - specified tolerance of eigenvalues.

  Niters=50^4;
  eigerr = zeros(Niters,1);
  errhist = zeros(Niters,1); 

  eigs2=copy(sort(eigs));

  H = copy(A);
  nn = size(H);
  n = nn[1];
  evals = zeros(neig,1);

# Going in Reverse order
  rend = n-neig+2;

#  eye = Matrix{Float64}(I,r,r);     # Identity matrix
  
  i = 0;
  j = 0;
  for r in n:-1:rend
    global subH
    global offdiag

    k=0;
    msize = r;          # Matrix size
    subH = H[1:r,1:r];
    offdiag = subH[r,r-1];
    eye = one(subH);  

    while offdiag > errtol
      k=k+1;
      j=j+1;

      if (ifrayleigh)
        shift=subH[r,r];
      end  

      subH = subH - shift*eye;
      ifmod = false;          # if modified GS
      ngs = 2;                # no of GS passes.

      Q,R = QRfac(subH,ifmod,ngs);
      subH = R*Q + eye*shift;
      approx_eigs = copy(sort(diag(subH)));
      
      eigerr[j] = maximum(abs.(approx_eigs - eigs2[1:r]));
      ii = argmax(abs.(approx_eigs - eigs2[1:r]));
      println([ii eigerr[j]])

      offdiag = abs(subH[r,r-1]);
     
      errhist[j]=geterr(subH);

      if nstep>0 && k==nstep
        offdiag=0.
      end
      
    end

    i=i+1;
    evals[i]=subH[r,r];  
    H[1:r,1:r] = subH;

  end

  return H,evals,errhist[1:j],eigerr[1:j]

end


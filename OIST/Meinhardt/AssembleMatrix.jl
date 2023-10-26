function AssembleMatrix(c0,cnv,ν,wlp,lx1,nel)

      dof = nel*(lx1-1) + 1;

      A   = zeros(Float64,dof,dof);

      for i in 1:nel
        j1 = (i-1)*(lx1-1) + 1;
        j2 = j1 + lx1 -1;
        subm = -c0.*cnv[:,:,i] + ν.*wlp[:,:,i];
        A[j1:j2,j1:j2] = A[j1:j2,j1:j2] + subm;
      end

      return A
end      

function SEM_QQT(glnum;Dtype=Float64)

      # QT          - Gather  (Direct Stiffness Summation)
      # Q           - Scatter (Global to local)
      # Dtype       - Precision (Float64/BigFloat)

      # glnum - Global numbering of points

      lx1,nel = size(glnum)

      ndof1 = nel*(lx1-1) + 1;
      ndof2 = nel*lx1

      Q     = spzeros(Dtype,ndof2,ndof1)
      QT    = spzeros(Dtype,ndof1,ndof2)

      # Built in a very rudimentary way.
      # Should be a more sophisticated way to build this matrix
      for i in 1:lx1, j in 1:nel
        # Corresponding global and local positions
        jg = glnum[i,j]
        jl = (j-1)*lx1 + i;

        QT[jg,jl] = QT[jg,jl] + 1.0
        Q[jl,jg]  = Q[jl,jg] + 1.0
      end

      return Q,QT
end

#---------------------------------------------------------------------- 

function SEM_Global_Num(xm1,ifperiodic::Bool;Dtype=Float64)

      lx1,nel = size(xm1)

      n = lx1*nel

      glnum = zeros(Int64,lx1,nel) 

      xtmp  = copy(xm1[:]);

      eps10 = eps(Dtype)   # Get machine precision of the variable type
      eps10 = 1000.0*eps10

      ii    = sortperm(xtmp)
      sort!(xtmp)

      gno   = 0
      y0    = Dtype(10000)
      x0    = xtmp[1] - y0
      diff  = 0.
      for k in 1:n
        ij = ii[k]
        i  = mod(ij,lx1)
        if (i==0)
          i = lx1
        end  
        j   = Int(floor((ij-i)/lx1)) + 1

        x1   = xtmp[k]
        diff = abs(x1 - x0)
        if diff < eps10
          glnum[i,j] = gno
        elseif i==lx1 && j == nel && ifperiodic
          glnum[i,j] = glnum[1,1]
        else
          gno        = gno + 1
          glnum[i,j] = gno
          x0         = x1
        end  
      end            
        
      return gno, glnum
    end      
#---------------------------------------------------------------------- 







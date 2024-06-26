function Sem_QQT(glnum,prec)

#     QT          - Gather  (Direct Stiffness Summation)
#     Q           - Scatter (Global to local)
#     prec        - Precision (Float64/BigFloat)

#     glnum - Global numbering of points

      lx1,nel = size(glnum)
      glnu    = unique(glnum)       # Unique Points

      ndof1 = length(glnu)
      ndof2 = nel*lx1

      Q     = spzeros(prec,ndof2,ndof1)
      QT    = spzeros(prec,ndof1,ndof2)

#     Built in a very rudimentary way.
#     Should be a more sophisticated way to build this matrix
      for i in 1:lx1, j in 1:nel
#       Corresponding global and local positions
        jg = glnum[i,j]
        jl = (j-1)*lx1 + i;

        QT[jg,jl] = QT[jg,jl] + 1.0
        Q[jl,jg]  = Q[jl,jg] + 1.0

      end

#      Q = transpose(QT)

      return Q,QT
end

#---------------------------------------------------------------------- 

function Sem_Global_Num(xm1,prec,ifp)

#     ifp - If end points are periodic?

      lx1,nel = size(xm1)

      n = lx1*nel

      glnum = zeros(Int64,lx1,nel) 

      xtmp  = xm1[:];
      if isinf(xtmp[1])
        xtmp[1] = prec(-999999999.999)
      end  

      if isinf(xtmp[n])
        xtmp[n] = prec(999999999.999)
      end  

      x1    = xtmp[1]
      eps10 = eps(prec)   # Get machine precision of the variable type
      eps10 = 1000.0*eps10

      ii    = sortperm(xtmp)
      sort!(xtmp)

      gno   = 0
      y0    = prec(10000)
      x0    = xtmp[1] - y0
      diff  = 0.
      for k in 1:n
            
        ij = ii[k]
        i  = mod(ij,lx1)
        if (i==0)
          i = lx1
        end  
        j  = Int(floor((ij-i)/lx1)) + 1

        x1 = xtmp[k]
        diff = abs(x1 - x0)
        if diff < eps10
          glnum[i,j] = gno
        elseif i==lx1 && j == nel && ifp
          glnum[i,j] = glnum[1,1]
        else
          gno = gno + 1
          glnum[i,j] = gno
          x0 = x1
        end  

      end            
        
      return gno, glnum
end      
#---------------------------------------------------------------------- 







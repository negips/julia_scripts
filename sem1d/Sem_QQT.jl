function Sem_QQT(glnum,prec)

#     QT          - Gather  (Direct Stiffness Summation)
#     Q           - Scatter (Global to local)
#     prec        - Precision (Float64/BigFloat)

#     glnum - Global numbering of points

      lx1,nel = size(glnum)

      ndof1 = nel*(lx1-1) + 1;
      ndof2 = nel*lx1

      Q     = zeros(prec,ndof2,ndof1)
      QT    = zeros(prec,ndof1,ndof2)

#     Built in a very rudimentary way.
#     Should be a more sophisticated way to build this matrix
      for i in 1:lx1, j in 1:nel
#       Corresponding global and local positions
        jg = glnum[i,j]
        jl = (j-1)*lx1 + i;

        QT[jg,jl] = QT[jg,jl] + 1.0

      end

      Q = transpose(QT)

      return Q,QT
end

#---------------------------------------------------------------------- 

function Sem_Global_Num(xm1,prec)

      lx1,nel = size(xm1)

      n = lx1*nel

      glnum = zeros(Int64,lx1,nel) 

      xtmp  = xm1[:];

      x1    = xtmp[1]
      eps10 = eps(prec)   # Get machine precision of the variable type
      eps10 = 1000.0*eps10

      ii    = sortperm(xtmp)
      sort!(xtmp)

      gno   = 0 
      if prec==BigFloat
        y0 = BigFloat(10000.0)
      else
        y0 = 10000.0
      end  
      x0    = xtmp[1] - y0
      diff  = 0.
      for k in 1:n
#        local x0::Float64
            
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
        else
          gno = gno + 1
          glnum[i,j] = gno
        x0 = x1
      end  


      end            
        
      return gno, glnum

    end      

#---------------------------------------------------------------------- 







function Sem_QQT(glnum)

#     QT - Gather  (Direct Stiffness Summation)
#     Q  - Scatter (Global to local)

#     glnum - Global numbering of points

      lx1,nel = size(glnum)

      ndof1 = nel*(lx1-1) + 1;
      ndof2 = nel*lx1

      Q     = zeros(Float64,ndof2,ndof1)
      QT    = zeros(Float64,ndof1,ndof2)

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

function Sem_Global_Num(xm1)

      lx1,nel = size(xm1)

      n = lx1*nel

      glnum = zeros(Int64,lx1,nel) 

      xtmp  = xm1[:];

      x1    = xtmp[1]
      eps10 = eps(typeof(x1))   # Get machine precision of the variable type
      eps10 = 100.0*eps10

      ii    = sortperm(xtmp)
      sort!(xtmp)

      gno   = 0
      x0    = xtmp[1] - 10000.
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







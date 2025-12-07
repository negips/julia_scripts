# Perform Explicit Shifted QR with exact eigenvalue shifts
function ExplicitShiftedQR(Hs,μ0,nμ,ngs)

  tol = 1.0e-16
  H   = copy(Hs)

  r,c = size(H)

  Q         = Matrix{eltype(H)}(I,r,c)    # Identity
  Qtmp      = Matrix{eltype(H)}(I,r,c)    # Identity

  R   = zeros(eltype(H),r,c)

  if nμ == 0
    nμ = 1
    μ  = [0.0]
  else
    μ  = μ0
  end  

  for i in 1:nμ
    k   = i
    Hj  = H - μ[k]*I
    Qj  = Matrix{eltype(Q)}(I,r,c)    # Identity
    q   = Hj[:,1]
    Qj[:,1] = q/norm(q)
    R[1,1]  = norm(q)
    for m in 2:r
      R[m,1]    = 0.
    end  
    for j in 2:c
      q  = Hj[:,j]
      h  = Qj[:,1:j-1]'*q
      q .= q .- Qj[:,1:j-1]*h
      if ngs>1
        for l = 1:ngs-1
          g  = Qj[:,1:j-1]'*q
          h  = h .+ g
          q .= q .- Qj[:,1:j-1]*g
        end
      elseif ngs<=0
        while true
          g   = Qj[:,1:j-1]'*q
          res = abs(g'*g)
          break
          if res<tol
            break
          end
          h = h .+ g
          q = q - Qj[:,1:j-1]*g
        end  
      end  
      Qj[:,j]     = q/norm(q)
      R[1:j-1,j]  = h
      R[j,j]      = norm(q)
      if j<r
        for m in j+1:r
          R[m,j]    = 0.
        end  
      end 
    end
   
    H        = R*Qj + μ[k]*I
    Q        = Qtmp*Qj
    Qtmp     = copy(Q)  

  end
 
  return H,Q 
end




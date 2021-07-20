# Perform Explicit Shifted QR with exact eigenvalue shifts
function ExplicitShiftedQR(Hs,μ0,nμ,ngs)

  tol = 1.0e-12
  H   = copy(Hs)

  r,c = size(H)

  Q   = Matrix{typeof(H[1,1])}(1.0I,r,c)    # Identity
  Qj  = Matrix{typeof(H[1,1])}(1.0I,r,c)    # Identity

  R   = copy(H)

  if nμ == 0
    nμ = 1
    μ  = [0.0]
  else
    μ  = μ0
  end  

  for i in 1:nμ
    Hj  = H - μ[i]*I
    Qj  = Matrix{typeof(Hj[1,1])}(1.0I,r,c)    # Identity
    q   = Hj[:,1]
    Qj[:,1] = q/norm(q)
    R[1,1]  = norm(q)
    R[2,1]  = 0.
    for j in 2:c
      q = Hj[:,j]
      h = Qj[:,1:j-1]'*q
      q .= q .- Qj[:,1:j-1]*h
      if ngs>1
        for k = 1:ngs-1
          g = Qj[:,1:j-1]'*q
          h = h .+ g
          q = q - Qj[:,1:j-1]*g
        end
      elseif ngs<=0
        while true
          g   = Qj[:,1:j-1]'*q
          res = abs(g'*g)
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
        R[j+1,j]    = 0.
      end 
    end
   
#    h = qr(H)
#    Qj = h.Q
#    R  = h.R
    H = R*Qj + μ[i]*I
    Q = Q*Qj

#    println(norm(Q'Q - I))
  end 
 
  return Q,H 
end




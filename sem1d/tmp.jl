      Hes = copy(Hold)
      VV  = Vold

      kk = LKryl
      k  = kk + 1
      r  = VV[:,k]*Hes[k,k-1]

      H  = H[1:kk,1:kk]

      F  = eigen(H)          # Uses Lapack routine (dgeev/zgeev)
      fr = real.(F.values)
      fr_sort_i = sortperm(fr,rev=false)   # Increasing order
      μ         = F.values[fr_sort_i[1:EKryl]]
      nμ        = length(μ)

      Q,Hs  = ExplicitShiftedQR(H,μ,nμ,0)
      v     = V[:,1:kk]*Q[:,Nev+1]        # Part of new residual vector
      βk    = Hs[Nev+1,Nev]               # e_k+1^T*H*e_k         # This in principle is zero
      σ     = Q[kk,Nev]                   # e_k+p^T*Q*e_k

      r2    = βk*v .+ σ*r                 # new residual vector
      β     = abs(sqrt(r2'*(Bg.*r2)))

      r2    = r2/β








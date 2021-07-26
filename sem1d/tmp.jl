      include("BulgeChase.jl")

      Hes = copy(Hold)
      VV  = Vold

      kk = LKryl
      k  = kk + 1
      r  = VV[:,k]*Hes[k,k-1]

      H  = Hes[1:kk,1:kk]

      F  = eigen(H)          # Uses Lapack routine (dgeev/zgeev)
      Fval = F.values
      fr = real.(F.values)
      fr_sort_i = sortperm(fr,rev=false)   # Increasing order
      μ         = F.values[fr_sort_i[1:EKryl]]
      nμ        = 50 #length(μ)

#      Hs,Q  = ExplicitShiftedQR(H,μ,nμ,2)
      Hs,Q  = FrancisSeq(H,μ,nμ)
      v     = V[:,1:kk]*Q[:,Nev+1]        # Part of new residual vector
      βk    = Hs[kk-nμ+1,kk-nμ]               # e_k+1^T*H*e_k         # This in principle is zero
      σ     = Q[kk,Nev]                   # e_k+p^T*Q*e_k

      r2    = βk*v .+ σ*r                 # new residual vector
      β     = abs(sqrt(r2'*(Bg.*r2)))

#      r2    = r2/β


      en          = zeros(ComplexF64,LKryl)
      en[LKryl]   = 2.0
      vn          = Vold[:,LKryl+1]*Hold[LKryl+1,LKryl]
      erM         = vn*en'

      erMQ        = erM*Q

      β_1         = Hs[kk-nμ-1,kk-nμ-2]
      β0          = Hs[kk-nμ,kk-nμ-1]
      β1          = Hs[kk-nμ+1,kk-nμ]
      β2          = Hs[kk-nμ+2,kk-nμ+1]

      θ           = [β_1 β0 β1 β2]

#     Reverse Step
      F2  = eigen(H)          # Uses Lapack routine (dgeev/zgeev)
      F2v = F2.values
      fr  = real.(F2.values)
      fr_sort_i = sortperm(fr,rev=true)   # Decreasing order
      μ2        = F2v[fr_sort_i[1:EKryl]]
      nμ2       = 50 #length(μ)

      Hs2,Q2  = RevFrancisSeq(H,μ2,nμ2)

      β_1         = Hs2[nμ2-1,nμ2-2]
      β0          = Hs2[nμ2,nμ2-1]
      β1          = Hs2[nμ2+1,nμ2]
      β2          = Hs2[nμ2+2,nμ2+1]

      θ2          = [β_1 β0 β1 β2]




# Arnoldi Implicit restart
function ArnIRst(V::Matrix,Hes::Matrix,B::Vector,k::Int,kmax::Int,Nev::Int,gs::Int)

#   V         - Krylov Vector
#   H         - Upper Hessenberg
#   B         - Weight (vector)
#   β         - previous/new residual norm
#   k         - Current Krylov size
#   kmax      - Maximum Krylov size
#   r         - Old/New residual vector for Arnoldi iteration
#   Nev       - Eigenvalues to retain
#   ngs       - No of Gram-Schmidt

    tol = 1.0e-12 

    EKryl = kmax - 1 - Nev 
    ifconv = false
    nkryl = k

    r = V[:,k]*Hes[k,k-1]

    rw,cl = size(V)

    U     = 0.0*V
    G     = 0.0*Hes

#   Perform implicit restart      
    if k == kmax

      kk = k-1

      H = Hes[1:kk,1:kk]

      F  = eigen(H)          # Uses Lapack routine (dgeev/zgeev)
      fr = real.(F.values)
      fr_sort_i = sortperm(fr,rev=false)   # Increasing order
      μ         = F.values[fr_sort_i[1:EKryl]]
      nμ        = length(μ)

      Q,Hs  = ExplicitShiftedQR(H,μ,nμ,ngs)
      v     = V[:,1:kk]*Q[:,Nev+1]        # Part of new residual vector
      βk    = Hs[Nev+1,Nev]               # e_k+1^T*H*e_k         # This in principle is zero
      σ     = Q[kk,Nev]                   # e_k+p^T*Q*e_k

      r    .= βk*v .+ σ*r                 # new residual vector
      β     = abs(sqrt(r'*(B.*r)))

      r     = r/β

      if abs(β)<tol
        ifconv = true
      end  

      U[:,1:Nev]        = V[:,1:kk]*Q[:,1:Nev]        # Updated Krylov space
      G[1:Nev,1:Nev]    = Hs[1:Nev,1:Nev]             # New Upper Hessenberg
      G[Nev+1,Nev]      = β
      U[:,Nev+1]        = r

      nkryl = Nev+1

    end     # k == kmax+1

    return U,G,nkryl,ifconv
end  






# Arnoldi Implicit restart
function ArnIRst!(V,H,B,k,kmax,v,ngs)
    global V,H,v

    EKryl = kmax - k

#   Perform implicit restart      
    if k == kmax+1
      global Hb
      kk = k-1
      Hb = H[1:kk,1:kk]
      F  = eigen(Hb)          # Uses Lapack routine (dgeev/zgeev)
      fr = real.(F.values)
      fr_sort_i = sortperm(fr,rev=false)   # Increasing order
      II = Matrix{Complex}(1.0I,kk,kk)
      QQ = deepcopy(II)
      for j in 1:1 #EKryl
        k = fr_sort_i[j]
        QRj = qr(Hb - F.values[k]*II)          # Expensive. Should replace by Bulge Chase
        Hb  = QRj.T*QRj.factors + F.values[k]*II
        QQ  = QQ*QRj.factors
      end 
      V = V*QQ
      β = Hb[Nev+1,Nev]       # e_k+1^T*Hb*e_k
      σ = QQ[kk,Nev]           # e_k+p^T*Q*e_k
      r = copy(v)
      v = V*QQ[:,Nev+1]        # V*Q*e_k+1
      r = β.*v .+ σ.*r        # Update residual

      H[1:kk,1:kk] = Hb

      β = sqrt(r'*(Bg.*r))
      v = r/β
      H[Nev+1,Nev] = β

      k            = Nev

      ifconv = true
    end     # k == kmax+1 



  return
end  


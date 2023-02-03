# BiOrtho Implicit restart
function BiOrthoIRst3!(V::Matrix,W::Matrix,Tv::Matrix,Tw::Matrix,fv::Vector,fw::Vector,k::Int,kmax::Int,Nev::Int,ngs::Int)

#   V         - Right Krylov Vector
#   W         - Left Krylov Vector
#   Tv        - Tridiagonal Matrix
#   Tw        - Tridiagonal Matrix (Slaved)
#   β         - previous/new residual norm
#   k         - Current Krylov size
#   kmax      - Maximum Krylov size
#   fv        - Old/New residual vector for Right subspace (scaled)
#   fw        - Old/New residual vector for Left subspace  (scaled)
#   Nev       - Eigenvalues to retain


    el        = eltype(V[1])
    tol       = 10000.0*eps(abs(V[1]))

    zro       = el(0.0) 

    ekryl = kmax - Nev 
    ifconv = false
    nkryl = k

    Vcopy = zro*copy(V)
    Wcopy = zro*copy(W)

    ekt     = zeros(el,1,k-1)
    ekt[k-1] = el(1.0)

#   Perform implicit restart 
    if k == kmax+1

      kk = k-1

      T = copy(Tv[1:kk,1:kk])

#     Perform QR operations 
#      μ,nμ   = GetLowerShifts(T,ekryl)                # Unwanted shifts
      μ,nμ   = GetLowerImaginaryShifts(T,ekryl)                # Unwanted shifts

#      Hsv,Qv = ExplicitShiftedQR(H,μ,nμ,ngs)
      Tp    = copy(T)
      Tpa   = copy(T')

#      v1,w1  = ImplicitLRSeq!(Tp,μ,nμ)
      ektv = copy(ekt)
      ektw = copy(ekt)
 
      Vcopy = V[:,1:kk]                     # Updated Right Krylov space
      Wcopy = W[:,1:kk]                     # Updated Right Krylov space

      for j in 1:nμ
        λ       = μ[j]
        wi,vi   = SimilarityTransformBulge!(Tp,λ,1)

        Vcopy .= Vcopy*vi
        Wcopy .= Wcopy*(wi')

        ektv  .= ektv*vi
        ektw  .= ektw*(wi')
 
        for i in 1:kk-2
          wi,vi  = SimilarityTransform!(Tp,i,kk)
          Vcopy .= Vcopy*vi
          Wcopy .= Wcopy*(wi')

          ektv  .= ektv*vi
          ektw  .= ektw*(wi')
        end

#       Adjoint Matrix        
        Tpa     = copy(Tp')
        λ       = μ[j]'
        wi,vi   = SimilarityTransformBulge!(Tpa,λ,1)

        Vcopy .= Vcopy*(wi')
        Wcopy .= Wcopy*vi

        ektv  .= ektv*(wi')
        ektw  .= ektw*vi
 
        for i in 1:kk-2
          wi,vi  = SimilarityTransform!(Tpa,i,kk)
          Vcopy .= Vcopy*(wi')
          Wcopy .= Wcopy*vi

          ektv  .= ektv*(wi')
          ektw  .= ektw*vi
        end

        Tp     = copy(Tpa')

      end   # j in 1:nμ 


#      Tpa    = copy(Tp')
#      μs     = adjoint.(μ)
#      v2,w2  = ImplicitLRSeq!(Tpa,μs,nμ)
#
#      Tp    .= Tpa'
#
#      # This is a row vector
#      ektv   = ekt*v1
#      ektvw  = ektv*(w2')
#
#      # This is a row vector
#      ektw   = ekt*(w1')
#      ektwv  = ektw*(v2)

      rv     = fv*ektv[Nev]
      rw     = fw*ektw[Nev] 

      βv     = Tp[Nev+1,Nev]                         # e_k+1^T*(Tp^T)*e_k         # This in principle is zero
      βw     = Tpa[Nev+1,Nev]                        # e_k+1^T*(Tp^T)*e_k         # This in principle is zero

#     Update Krylov spaces      
#      Vcopy             = V[:,1:kk]*(v1*w2')                     # Updated Right Krylov space
      V                .= zro.*V
      V[:,1:Nev]       .= Vcopy[:,1:Nev]
      rV                = Vcopy[:,Nev+1]

#      Vcopy             = W[:,1:kk]*(w1'*v2)                     # Updated Left Krylov space
      W                .= zro.*W
      W[:,1:Nev]       .= Wcopy[:,1:Nev]
      rW                = Wcopy[:,Nev+1]

#     Update Triangular Matrices
      Tv               .= zro.*Tv     
      Tv[1:Nev,1:Nev]   = Tp[1:Nev,1:Nev]                   # New Tridiagonal Matrix

      Tw               .= zro.*Tw
      Tw[1:Nev,1:Nev]   = Tpa[1:Nev,1:Nev]                  # New Tridiagonal Matrix

#      @printf "βv After ImplicitLR: %12e, %12eim\n" real(βv) imag(βv) 

      fv    .= rv .+ βv*rV                     # new right residual vector
      fw    .= rw .+ βw*rW                     # new left  residual vector

      θ       = fw'*fv                         # <̂w,̂v>
      θa      = abs(θ)

      vn      = norm(fv)
      wn      = norm(fw)

      if vn<tol
        println("Possible Right subspace convergence ||v|| = $vn")
        ifconv = true
      end  

      if wn<tol
        println("Possible Left subspace convergence ||w|| = $wn")
        ifconv = true
      end  

      if abs(θ)<tol
        println("Possible Orthogonal subspaces: |<w,v>| = $θa, ||v|| = $vn, ||w|| = $wn ")
        ifconv = true
      end  

      if (ifconv)
        V[:,Nev+1]      = fv
        W[:,Nev+1]      = fw
        
        nkryl = Nev+1

        return nkryl,ifconv
      end  

      δv,δw = BiorthoScale_vw!(fv,fw)
#      println("θ  = $θ")
#      println("δv = $δv")
#      println("δw = $δw")

#     Update Krylov spaces 
      V[:,Nev+1]        = fv
      W[:,Nev+1]        = fw

#     Update Hessenberg Matrices      
      Tv[Nev+1,Nev]     = δv
      Tw[Nev+1,Nev]     = δw

      nkryl = Nev+1

    else
      ifconv            = false
      nkryl             = k

    end     # k == kmax+1

    return nkryl,ifconv
end  
#----------------------------------------------------------------------

# BiOrtho Implicit restart
function BiOrthoIRst2!(V::Matrix,W::Matrix,Tv::Matrix,Tw::Matrix,fv::Vector,fw::Vector,k::Int,kmax::Int,Nev::Int,ngs::Int)

#   V         - Right Krylov Vector
#   W         - Left Krylov Vector
#   Tv        - Tridiagonal Matrix
#   Tw        - Tridiagonal Matrix (Slaved)
#   β         - previous/new residual norm
#   k         - Current Krylov size
#   kmax      - Maximum Krylov size
#   fv        - Old/New residual vector for Right subspace (scaled)
#   fw        - Old/New residual vector for Left subspace  (scaled)
#   Nev       - Eigenvalues to retain


    el        = eltype(V[1])
    tol       = 10000.0*eps(abs(V[1]))

    zro       = el(0.0) 

    ekryl = kmax - Nev 
    ifconv = false
    nkryl = k

    Vcopy = zro*copy(V)
    Wcopy = zro*copy(W)

    ekt     = zeros(el,1,k-1)
    ekt[k-1] = el(1.0)

#   Perform implicit restart 
    if k == kmax+1

      kk = k-1

      T = copy(Tv[1:kk,1:kk])

#     Perform QR operations 
      μ,nμ   = GetLowerShifts(T,ekryl)                # Unwanted shifts
#      Hsv,Qv = ExplicitShiftedQR(H,μ,nμ,ngs)
      μ,nμ   = GetLowerImaginaryShifts(T,ekryl)                # Unwanted shifts

      Tp    = copy(T)      
      v1,w1  = ImplicitLRSeq!(Tp,μ,nμ)

      Tpa    = copy(Tp')
      μs     = adjoint.(μ)
      v2,w2  = ImplicitLRSeq!(Tpa,μs,nμ)

      Tp    .= Tpa'

      # This is a row vector
      ektv   = ekt*v1
      ektvw  = ektv*(w2')

      # This is a row vector
      ektw   = ekt*(w1')
      ektwv  = ektw*(v2)

#     If we build the full residual matrix      
#      rMv    = fv*ektv
      rv     = fv*ektvw[Nev]

#     If we build the full residual matrix      
#      rMw    = fw*ektw
      rw     = fw*ektwv[Nev] 

      βv     = Tp[Nev+1,Nev]                         # e_k+1^T*(Tp^T)*e_k         # This in principle is zero
      βw     = Tpa[Nev+1,Nev]                        # e_k+1^T*(Tp^T)*e_k         # This in principle is zero

      red    = (v2'*w1)*(v1*w2')

      redn1   = norm(w1*v1 - I)
      redn2   = norm(v2'*w2' - I)
      rednorm = norm(red - I)

#      if redn2>1.0e-8
#        display(Tp)
#      end  

#     Update Krylov spaces      
      Vcopy             = V[:,1:kk]*(v1*w2')                     # Updated Right Krylov space
      V                .= zro.*V
      V[:,1:Nev]       .= Vcopy[:,1:Nev]
      rV                = Vcopy[:,Nev+1]

      invm              = inv(v1*w2')

      Wcopy             = W[:,1:kk]*(w1'*v2)                     # Updated Left Krylov space
#      Wcopy             = W[:,1:kk]*(invm')                       # Updated Left Krylov space
      W                .= zro.*W
      W[:,1:Nev]       .= Wcopy[:,1:Nev]
      rW                = Wcopy[:,Nev+1]

#     Update Triangular Matrices
      Tv               .= zro.*Tv     
      Tv[1:Nev,1:Nev]   = Tp[1:Nev,1:Nev]                   # New Tridiagonal Matrix

      Tw               .= zro.*Tw
      Tw[1:Nev,1:Nev]   = Tpa[1:Nev,1:Nev]                  # New Tridiagonal Matrix

#      @printf "βv After ImplicitLR: %12e, %12eim\n" real(βv) imag(βv) 

      fv    .= rv .+ βv*rV                     # new right residual vector
      fw    .= rw .+ βw*rW                     # new left  residual vector

      fvo    = norm(W'*fv)
      fwo    = norm(V'*fw)

      βvn    = abs(βv)
      βwn    = abs(βw)

#      println("LU norm = $redn1, $redn2 $rednorm $fvo $fwo $βvn $βwn")

      θ       = fw'*fv                                # <̂w,̂v>
      θa      = abs(θ)

      vn      = norm(fv)
      wn      = norm(fw)

      if vn<tol
        println("Possible Right subspace convergence ||v|| = $vn")
        ifconv = true
      end  

      if wn<tol
        println("Possible Left subspace convergence ||w|| = $wn")
        ifconv = true
      end  

      if abs(θ)<tol
        println("Possible Orthogonal subspaces: |<w,v>| = $θa, ||v|| = $vn, ||w|| = $wn ")
        ifconv = true
      end  

      if (ifconv)
        V[:,Nev+1]      = fv
        W[:,Nev+1]      = fw
        
        nkryl = Nev+1

        return nkryl,ifconv
      end  

      δv,δw = BiorthoScale_vw!(fv,fw)
#      println("θ  = $θ")
#      println("δv = $δv")
#      println("δw = $δw")

#     Update Krylov spaces 
      V[:,Nev+1]        = fv
      W[:,Nev+1]        = fw

#     Update Hessenberg Matrices      
      Tv[Nev+1,Nev]     = δv
      Tw[Nev+1,Nev]     = δw

      nkryl = Nev+1

    else
      ifconv            = false
      nkryl             = k

    end     # k == kmax+1

    return nkryl,ifconv
end  
#----------------------------------------------------------------------
function BiOrthoIRst!(V::Matrix,W::Matrix,Hv::Matrix,Hw::Matrix,fv::Vector,fw::Vector,k::Int,kmax::Int,Nev::Int,ngs::Int)

#   V         - Right Krylov Vector
#   W         - Left Krylov Vector
#   Hv        - Right Upper Hessenberg
#   Hw        - Left Upper Hessenberg
#   β         - previous/new residual norm
#   k         - Current Krylov size
#   kmax      - Maximum Krylov size
#   fv        - Old/New residual vector for Right subspace
#   fw        - Old/New residual vector for Left subspace
#   Nev       - Eigenvalues to retain


    el        = eltype(V[1])
    tol       = eps(abs(V[1]))

    zro       = el(0.0) 

    ekryl = kmax - Nev 
    ifconv = false
    nkryl = k


    rv = V[:,k]*Hv[k,k-1]
    rw = W[:,k]*Hw[k,k-1]

    Vcopy = zro*V
    Wcopy = zro*W
    Hcopy = zro*Hv

#   Perform implicit restart      
    if k == kmax+1

      kk = k-1

      H = copy(Hv[1:kk,1:kk])

#     Perform QR operations      
      μ,nμ   = GetLowerShifts(H,ekryl)
      Hsv,Qv  = ExplicitShiftedQR(H,μ,nμ,ngs)

#      Hs,Q  = FrancisSeq(H,b,μ,nμ)
#      Hs,Q  = FrancisSeqExact(H,b,μ,nμ)

      v     = V[:,1:kk]*Qv[:,Nev+1]        # Part of new residual vector
      βv    = Hsv[Nev+1,Nev]               # e_k+1^T*H*e_k         # This in principle is zero

      @printf "βv After ImplicitQR: %12e, %12eim\n" real(βv) imag(βv) 

      H = copy(Hw[1:kk,1:kk])

#     Perform QR operations      
      μ,nμ    = GetLowerShifts(H,ekryl)
      Hsw,Qw  = ExplicitShiftedQR(H,μ,nμ,ngs)

      w     = W[:,1:kk]*Qw[:,Nev+1]       # Part of new residual vector
      βw    = Hsw[Nev+1,Nev]               # e_k+1^T*H*e_k         # This in principle is zero

      @printf "βw After ImplicitQR: %12e, %12eim\n" real(βw) imag(βw) 

#      hdiff = norm(H) - norm(Hs);
#      if (abs(hdiff)>1.0e-12)
#        @printf "Possible Forward instability: HDiff: %8e\n" hdiff
#      end  

      σv     = Qv[kk,Nev]                   # e_k+p^T*Q*e_k
      fv    .= βv*v .+ σv*rv                # new residual vector

      σw     = Qw[kk,Nev]                   # e_k+p^T*Q*e_k
      fw    .= βw*w .+ σw*rw                # new residual vector

      v1      = copy(fv)
      w1      = copy(fw)

      θ       = fw'*fv                   # <̂w,̂v>
      θa      = abs(θ) 
      δ       = sqrt(θa)
      β       = (θ/δ)'

      fv     .= fv./δ
      fw     .= fw./β

      δv      = (fv'*v1)/(fv'*fv)
      δw      = (fw'*w1)/(fw'*fw)

      if abs(θ)<tol
        println("Possible Orthogonal subspaces: <w,v> = $θa")
        ifconv = true
      end  

#     Update Krylov spaces      
      Vcopy[:,1:Nev]    = V[:,1:kk]*Qv[:,1:Nev]        # Updated Right Krylov space
      V                .= zro.*V
      V[:,1:Nev]        = Vcopy[:,1:Nev]
      V[:,Nev+1]        = fv

      Vcopy[:,1:Nev]    = W[:,1:kk]*Qw[:,1:Nev]       # Updated Left Krylov space
      W                .= zro.*W
      W[:,1:Nev]        = Vcopy[:,1:Nev]
      W[:,Nev+1]        = fw

#     Update Hessenberg Matrices      
      Hcopy[1:Nev,1:Nev]  = Hsv[1:Nev,1:Nev]           # New Upper Hessenberg
      Hv                 .= zro*Hv
      Hv[1:Nev,1:Nev]    .= Hcopy[1:Nev,1:Nev]
      Hv[Nev+1,Nev]       = δv

      Hcopy[1:Nev,1:Nev]  = Hsw[1:Nev,1:Nev]           # New Upper Hessenberg
      Hw                 .= zro*Hw
      Hw[1:Nev,1:Nev]    .= Hsw[1:Nev,1:Nev]
      Hw[Nev+1,Nev]       = δw

      nkryl = Nev+1

    else
      ifconv      = false
      nkryl       = k

    end     # k == kmax+1

    return nkryl,ifconv
end  
#----------------------------------------------------------------------

function GetUpperShifts(H::Matrix,Nev::Int)

      r,c = size(H)
      f  = eigvals(H)          # Uses Lapack routine (dgeev/zgeev)
      fr = real.(f)
      fr_sort_i = sortperm(fr,rev=true)   # Decreasing order
      μ         = f[fr_sort_i[1:Nev]]
      nμ        = length(μ)

      return μ,nμ
end

#----------------------------------------------------------------------
function GetLowerShifts(H::Matrix,EKryl::Int)

      r,c = size(H)
      f  = eigvals(H)          # Uses Lapack routine (dgeev/zgeev)
      fr = real.(f)
      fr_sort_i = sortperm(fr,rev=false)   # Increasing order
      μ0        = f[fr_sort_i[1:EKryl]]
      nμ        = length(μ0)
#      μ         = μ0[nμ:-1:1]
      μ         = μ0      

      return μ,nμ
end

#----------------------------------------------------------------------
function GetLowerImaginaryShifts(H::Matrix,EKryl::Int)

      r,c = size(H)
      f  = eigvals(H)          # Uses Lapack routine (dgeev/zgeev)
      fr = real.(f)
      fi = imag.(f)
      fi_sort_i = sortperm(fi,rev=false)   # Increasing order
      μ0        = f[fi_sort_i[1:EKryl]]
      nμ        = length(μ0)
#      μ         = μ0[nμ:-1:1]
      μ         = μ0      

      return μ,nμ
end

#----------------------------------------------------------------------

#----------------------------------------------------------------------










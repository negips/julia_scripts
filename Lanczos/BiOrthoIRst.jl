# BiOrtho Implicit restart
function BiOrthoIRst2!(V::Matrix,W::Matrix,Hv::Matrix,Hw::Matrix,fv::Vector,fw::Vector,k::Int,kmax::Int,Nev::Int,ngs::Int)

#   V         - Right Krylov Vector
#   W         - Left Krylov Vector
#   Hv        - Right Upper Hessenberg
#   Hw        - Left Upper Hessenberg
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

    ek      = zeros(el,k-1)
    ek[k-1] = el(1.0)

#   Perform implicit restart 
    if k == kmax+1

      kk = k-1

      H = copy(Hv[1:kk,1:kk])

#     Perform QR operations 
      μ,nμ   = GetLowerShifts(H,ekryl)                # Unwanted shifts
      Hsv,Qv = ExplicitShiftedQR(H,μ,nμ,ngs)

      ektQ   = transpose(Qv[kk,:]);                   # Making sure this is a row vector
      rMv    = fv*ektQ
      rMw    = fw*ektQ

      βv     = Hsv[Nev+1,Nev]                         # e_k+1^T*H*e_k         # This in principle is zero

#     Update Krylov spaces      
      Vcopy[:,1:Nev]    = V[:,1:kk]*Qv[:,1:Nev]       # Updated Right Krylov space
      V                .= zro.*V
      V[:,1:Nev]        = Vcopy[:,1:Nev]

      Vcopy[:,1:Nev]    = W[:,1:kk]*Qv[:,1:Nev]       # Updated Left Krylov space
      W                .= zro.*W
      W[:,1:Nev]        = Vcopy[:,1:Nev]

#     Update Hessenberg Matrices      
      H[1:Nev,1:Nev]    = Hsv[1:Nev,1:Nev]            # New Upper Hessenberg
      Hv               .= zro*Hv
      Hv[1:Nev,1:Nev]  .= H[1:Nev,1:Nev]

      Hw               .= zro*Hw
      Hw[1:Nev,1:Nev]  .= (Hsv[1:Nev,1:Nev])'

#      @printf "βv After ImplicitQR: %12e, %12eim\n" real(βv) imag(βv) 

      fv    .= rMv[:,Nev]                             # new right residual vector
      fw    .= rMw[:,Nev]                             # new left  residual vector

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

#     Update Krylov spaces 
      V[:,Nev+1]        = fv
      W[:,Nev+1]        = fw

#     Update Hessenberg Matrices      
      Hv[Nev+1,Nev]     = δv
      Hw[Nev,Nev+1]     = δw

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

function LutzCustomShifts(H::Matrix,Nev::Int)

#     Analytical Eigenvalues
      ω1 = find_zero(airyai,(-3.0,-0.0))
      ω2 = find_zero(airyai,(-5.0,-3.0))
      ω3 = find_zero(airyai,(-6.0,-5.0))
      ω4 = find_zero(airyai,(-7.0,-6.0))
      ω5 = find_zero(airyai,(-8.0,-7.0))
      ω6 = find_zero(airyai,(-9.5,-8.0))
      ω7 = find_zero(airyai,(-10.5,-9.5))
      ω8 = find_zero(airyai,(-11.8,-10.5))
      ω9 = find_zero(airyai,(-12.0,-11.8))
      ω10 = find_zero(airyai,(-12.9,-12.0))
      ω11 = find_zero(airyai,(-13.8,-12.9))
      ω12 = find_zero(airyai,(-14.8,-13.8))
      ω13 = find_zero(airyai,(-15.8,-14.8))
      ω14 = find_zero(airyai,(-16.8,-15.8))
      ω15 = find_zero(airyai,(-17.5,-16.8))
      
      ω  = [ω1, ω2, ω3, ω4, ω5, ω6, ω7, ω8, ω9, ω10, ω11, ω12, ω13, ω14, ω15]
      U  = 6.0
      γ  = 1.0 - im*1.0
      
      Ω  = im*(U*U/8.0 .- U*U/(4.0*γ) .+ γ^(1.0/3.0)*(U^(4.0/3.0))/(160.0^(2.0/3.0))*ω)

      i         = 1
      μ         = Ω[i:i+Nev-1]
      nμ        = length(μ)

      return μ,nμ
end

#----------------------------------------------------------------------










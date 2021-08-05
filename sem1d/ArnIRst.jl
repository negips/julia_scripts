# Arnoldi Implicit restart
function ArnIRst(V::Matrix,Hes::Matrix,B::Union{Vector,Matrix},k::Int,kmax::Int,Nev::Int,ngs::Int)

#   V         - Krylov Vector
#   H         - Upper Hessenberg
#   B         - Weight (vector)
#   β         - previous/new residual norm
#   k         - Current Krylov size
#   kmax      - Maximum Krylov size
#   r         - Old/New residual vector for Arnoldi iteration
#   Nev       - Eigenvalues to retain
#   ngs       - No of Gram-Schmidt

    revFrancis = false 

    tol = 1.0e-16

    EKryl = kmax - 1 - Nev 
    ifconv = false
    nkryl = k

    r1 = V[:,k]*Hes[k,k-1]

    rw,cl = size(V)

    U     = 0.0*V
    G     = 0.0*Hes

#   Perform implicit restart      
    if k == kmax

      kk = k-1

      H = deepcopy(Hes[1:kk,1:kk])

      if (revFrancis)

#        μ,nμ  = ArnGetCustomShifts(H,Nev)
        μ,nμ  = ArnGetUpperShifts(H,Nev)

        Hs,Q  = RevFrancisSeq(H,μ,nμ)     
        v     = V[:,1:kk]*Q[:,Nev+1]        # Part of new residual vector
        βk    = Hs[Nev+1,Nev]               # e_k+1^T*H*e_k         # This in principle is zero
        println("βk After ImplicitQR: $βk")
        σ     = Q[kk,Nev]                   # e_k+p^T*Q*e_k

        r     = βk*v .+ σ*r1                # new residual vector
        if ndims(B)>1
          β     = abs(sqrt(r'*(B*r)))
        else
          β     = abs(sqrt(r'*(B.*r)))
        end  

        r     = r/β

      else

        μ,nμ  = ArnGetLowerShifts(H,EKryl)

        Hs,Q  = ExplicitShiftedQR(H,μ,nμ,ngs)
#        Hs,Q  = FrancisSeq(H,μ,nμ)
#        Hs,Q  = FrancisSeqExact(H,μ,nμ)     
        v     = V[:,1:kk]*Q[:,Nev+1]        # Part of new residual vector
        βk    = Hs[Nev+1,Nev]               # e_k+1^T*H*e_k         # This in principle is zero
        println("βk After ImplicitQR: $βk")

#        ResM  = r*Q[kk,:]
#        r     = ResM[:,Nev] 

        σ     = Q[kk,Nev]                   # e_k+p^T*Q*e_k
        r     = βk*v .+ σ*r1                 # new residual vector

        if ndims(B)>1
          β     = abs(sqrt(r'*(B*r)))
        else
          β     = abs(sqrt(r'*(B.*r)))
        end  

        r     = r/β
      end        

      if abs(β)<tol
        ifconv = true
      end  

      U[:,1:Nev]        = V[:,1:kk]*Q[:,1:Nev]        # Updated Krylov space
      G[1:Nev,1:Nev]    = Hs[1:Nev,1:Nev]             # New Upper Hessenberg
      G[Nev+1,Nev]      = β
      U[:,Nev+1]        = r

      nkryl = Nev+1

    else
      U           = V
      G           = Hes
      ifconv      = false
      nkryl       = k

    end     # k == kmax+1

    return U,G,nkryl,ifconv
end  
#----------------------------------------------------------------------

function ArnGetUpperShifts(H::Matrix,Nev::Int)

      r,c = size(H)
      f  = eigvals(H)          # Uses Lapack routine (dgeev/zgeev)
      fr = real.(f)
      fr_sort_i = sortperm(fr,rev=true)   # Decreasing order
      μ         = f[fr_sort_i[1:Nev]]
      nμ        = length(μ)

      return μ,nμ
end

#----------------------------------------------------------------------
function ArnGetLowerShifts(H::Matrix,EKryl::Int)

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

function ArnGetCustomShifts(H::Matrix,Nev::Int)

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










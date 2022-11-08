function IRLanczos!(Vin::Matrix,Win::Matrix,Tm::Matrix,AHw::Vector,Av::Vector,γ::Vector,k::Int,kmax::Int,Mi::Int,Nev::Int)

#   Vin     - Right Krylov Matrix
#   Win     - Left Krylov Matrix
#   Hes     - Hessenberg Matrix
#   B       - Weight Vector/Matrix
#   w       - New left vector
#   v       - New right vector
#   γ       - α,β,δ
#   k       - Current Krylov Size
#   kmax    - Max Krylov size
#   Mi      - Major iteration Number
#   Nev     - No of Requested Eigenvalues/Vectors

    V = Vin
    H = Hes

#   Update Arnoldi Vector
    LanczosUpd!(Au,AHw,V,W,γ,k)
    k             = k +1
    α   = γ[1]
    β   = γ[2]
    δ   = γ[3]

    j = k
    if j==LKryl
      Tj[j,j]      = α
    else
      Tj[j,j]      = α
      Tj[j+1,j]    = δ
      Tj[j,j+1]    = β
    end  

    j = j+1
    if (j<=kmax+1)
      Q[:,nkryl]     = Au
      P[:,nkryl]     = AHw
      mi             = Mi
    end 

#   Perform implicit restart      
    if j == kmax+1
      U,G,k2,ifconv = ArnIRst(V,H,b,B,k,kmax+b,Nev,ngs)
      V = U
      H = G
      
      v = V[:,k2]
      β = abs(H[Nev+1,Nev])
      @printf "Major Iteration: %3i; β: %8e\n" Mi β

      nkryl = k2
      mi    = mi + 1
    end  

    return nkryl,β,mi 

end  

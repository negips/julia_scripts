function IRAM!(Vin::Matrix,Hes::Matrix,B::Union{Vector,Matrix},v::Vector,k::Int,kmax::Int,Mi::Int,Nev,ngs::Int)

#   Vin     - Krylov Matrix
#   Hes     - Hessenberg Matrix
#   B       - Weight Vector/Matrix
#   v       - New Arnoldi vector
#   k       - Current Krylov Size
#   kmax    - Max Krylov size
#   Mi      - Major iteration Number
#   Nev     - No of Requested Eigenvalues/Vectors
#   ngs     - No of Gram-Schmidt orthogonalizations

    V = Vin
    H = Hes

    b = 1         # Block size
#   Update Arnoldi Vector
    h,β,r         = ArnUpd(V,b,B,v,k,ngs)
    k             = k +1
    if (k<=kmax+b)
      V[:,k]      = r
      nkryl       = k
      mi          = Mi
    end        
    nkryl         = k

    if (k>b)
      kb          = k-b
      k1          = k-1
      H[1:k1,kb]  = h
      H[k,kb]     = β
    end  

#   Perform implicit restart      
    if k == kmax+b
      U,G,k2,ifconv = ArnIRst(V,H,b,B,k,kmax+b,Nev,ngs)
      V = U
      H = G
      
      v = V[:,k2]
      β = abs(H[Nev+1,Nev])
      @printf "Major Iteration: %3i; β: %8e\n" Mi β

      nkryl = k2
      mi    = mi + 1
    end  

    return V,H,nkryl,β,mi 

end  

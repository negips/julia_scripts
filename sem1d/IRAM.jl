function IRAM!(Vin::Matrix,Hes::Matrix,Bg::Vector,v::Vector,k::Int,kmax::Int,Mi::Int,Nev,ngs::Int)

#   Vin     - Krylov Matrix
#   Hes     - Hessenberg Matrix
#   Bg      - Weight Vector
#   v       - New Arnoldi vector
#   k       - Current Krylov Size
#   kmax    - Max Krylov size
#   Mi      - Major iteration Number
#   Nev     - No of Requested Eigenvalues/Vectors
#   ngs     - No of Gram-Schmidt orthogonalizations
#

    V = Vin
    H = Hes

#   Update Arnoldi Vector
    if k == 0
      h,β,r = ArnUpd(V,Bg,v,k,ngs)
      V[:,1]        = r
      nkryl         = 1
      mi            = Mi  
    else
      h,β,r = ArnUpd(V,Bg,v,k,ngs)
      H[1:k,k]      = h
      H[k+1,k]      = β
      k             = k + 1
      V[:,k]        = r
      v             = r
      nkryl         = k
      mi            = Mi
    end  

#   Perform implicit restart      
    if k == kmax+1
      Hold = H
      Vold = V
      U,G,k2,ifconv = ArnIRst(V,H,Bg,k,kmax+1,Nev,ngs)
      V = U
      H = G
      
      v = V[:,k2]
      β = abs(H[Nev+1,Nev])
      println("Outer Iteration: $Mi; β=$β")

      nkryl = k2
      mi    = mi + 1
    end  


    return V,H,nkryl,β,mi 

end  

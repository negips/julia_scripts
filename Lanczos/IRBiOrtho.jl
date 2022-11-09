function IRBiOrtho!(Vin::Matrix,Win::Matrix,Hv::Matrix,Hw::Matrix,Av::Vector,AHw::Vector,k::Int,kmax::Int,Mi::Int,Nev::Int,ngs::Int)

# V         - Right Krylov Space
# W         - Left Krylov Space
# Av        - New vector Av     for k>0
# AHw       - New vector (A^H)w for k>0 
# k         - Current Krylov size
# kmax      - Maximum Krylov size
# Mi        - Major iteration Number
# Nev       - No of Requested Eigenvalues/Vectors
# ngs       - No of Gram-Schmidt orthogonalizations

    el  = eltype(Av)
    zro = el(0.0)

    γv .= zeros(el,kmax+1)
    γw .= zeros(el,kmax+1)

#   Update Biorthogonal Vectors
    BiOrthoUpd!(Av,AHw,V,W,γv,γw,k)

    Hv[:,k] = γv
    Hw[:,k] = γw
#    Hw[k,:] = γw


    V[:,k+1]  = Av
    W[:,k+1]  = AHw
    k2        = k+1

    mi        = Mi

#   Perform implicit restart      
    if k2 == kmax+1
      k2,ifconv = BiOrthoIRst2!(V,W,Hv,Hw,Av,AHw,k2,kmax,Nev,ngs)
      
      v = V[:,k2]
      β = abs(Hv[Nev+1,Nev])
      @printf "Major Iteration: %3i; β: %8e\n" Mi β

      mi    = mi + 1
    end  

    nkryl = k2  

    return nkryl,mi 

end  

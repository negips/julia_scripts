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

    ifconv = false
    el  = eltype(Av)
    zro = el(0.0)

    γv = zeros(el,kmax+1)
    γw = zeros(el,kmax+1)

#   Update Biorthogonal Vectors
    BiOrthoUpd!(Av,AHw,V,W,γv,γw,k,ngs)

    Hv[:,k] = γv
    Hw[:,k] = γw


    V[:,k+1]  = Av
    W[:,k+1]  = AHw
    k2        = k+1

    mi        = Mi

#   Perform implicit restart      
    if k2 == kmax+1

      Av    .= Hv[k2,k]*Av
      AHw   .= Hw[k2,k]*AHw
      k2,ifconv = BiOrthoIRst2!(V,W,Hv,Hw,Av,AHw,k2,kmax,Nev,ngs)
      
      β = abs(Hv[Nev+1,Nev])
      @printf "Major Iteration: %3i; β: %8e\n" Mi β

      mi    = mi + 1
    end  

    nkryl = k2  

    return nkryl,mi,ifconv 

end
#---------------------------------------------------------------------- 

function BiorthoReortho!(V::Matrix,W::Matrix,j::Int,ngs::Int)

# V         - Right Krylov Space
# W         - Left Krylov Space
# γv        - Oblique Projections onto V along W
# γw        - Oblique Projections onto W along V
# j         - Current Krylov size
# jmax      - Maximum Krylov size

  el  = eltype(u)
  zro = el(0.0)

  γr = zeros(el,j)
  γl = zeros(el,j)

# New Lanczos Vectors
  for l in 1:j
    for g = 1:ngs
      for k in 1:l-1
        γr[k]  = W[:,k]'*V[:,l]
      end
  
      for k in 1:l-1
        V[:,l] .= V[:,l] .- γr[k]*V[:,k]      # ̂v =  A*v_j+1   - β_j*v
      end 
    end     # g =1:ngs 
    
    v1      = copy(V[:,l])
    w1      = copy(AHw)

    θ       = W[:,l]'*V[:,l]            # <̂w,̂v>
    δ       = sqrt(abs(θ))
    β       = (θ/δ)'

    V[:,l] .= Av./δ
    W[:,l] .= AHw./β

  end  

end



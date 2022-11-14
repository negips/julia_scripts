# Lanczos Update
function BiOrthoUpd!(Av::Vector,AHw::Vector,V::Matrix,W::Matrix,γv::Vector,γw::Vector,j::Int,ngs::Int)

# Yousef Saad - Numerical Methods for large eigenvalue problems (2011) 2nd Edition
# Algorithm 6.6: The non-Hermitian Lanczos Algorithm

# V         - Right Krylov Space
# W         - Left Krylov Space
# Av        - New vector Av     for k>0
# AHw       - New vector (A^H)w for k>0 
# γv        - Oblique Projections onto V along W
# γw        - Oblique Projections onto W along V
# j         - Current Krylov size
# jmax      - Maximum Krylov size
# ngs       - No of oblique Gram-Schmid procedures

  el  = eltype(Av[1])
  zro = el(0.0)

  γv .= zro.*γv
  γw .= zro.*γw

  γr = zro.*γv
  γl = zro.*γw

# New Lanczos Vectors
  if j > 0

    for g = 1:ngs
      for k in 1:j
        γr[k]  = W[:,k]'*Av
        γl[k]  = V[:,k]'*AHw
      end
  
      for k in 1:j
        Av    .= Av  .- γr[k]*V[:,k]      # ̂v =  A*v_k+1   - γv_k*v
        AHw   .= AHw .- γl[k]*W[:,k]      # ̂w = (A')w_k+1  - γw_k*w
      end 
      γv    .= γv .+ γr
      γw    .= γw .+ γl
    end     # g =1:ngs 

#   Normalized s.t. <w,v> = 1
    αv,αw = BiorthoScale_vw!(Av,AHw)
    γv[j+1] = αv
    γw[j+1] = αw

  else

#   Normalized s.t. <w,v> = 1
    αv,αw = BiorthoScale_vw!(Av,AHw)
    γv[j+1] = αv
    γw[j+1] = αw
   
  end

end

#----------------------------------------------------------------------
# Lanczos Update with Orthogonal W
function BiOrthoUpd2!(Av::Vector,AHw::Vector,V::Matrix,W::Matrix,γv::Vector,γw::Vector,j::Int,ngs::Int)

# Yousef Saad - Numerical Methods for large eigenvalue problems (2011) 2nd Edition
# Algorithm 6.6: The non-Hermitian Lanczos Algorithm

# V         - Right Krylov Space
# W         - Left Krylov Space
# Av        - New vector Av     for k>0
# AHw       - New vector (A^H)w for k>0 
# γv        - Oblique Projections onto V along W
# γw        - Oblique Projections onto W along V
# j         - Current Krylov size
# jmax      - Maximum Krylov size
# ngs       - No of oblique Gram-Schmid procedures

  el  = eltype(Av[1])
  zro = el(0.0)

  γv .= zro.*γv
  γw .= zro.*γw

  γr = zro.*γv
  γl = zro.*γw

# New Lanczos Vectors
  if j > 0
    for g = 1:ngs
      for k in 1:j
        γr[k]  = W[:,k]'*Av
        γl[k]  = W[:,k]'*AHw
      end
  
      for k in 1:j
        Av    .= Av  .- γr[k]*V[:,k]      # ̂v =  A*v_k+1   - γv_k*v
        AHw   .= AHw .- γl[k]*W[:,k]      # ̂w = (A')w_k+1  - γw_k*w
      end 
      γv    .= γv .+ γr
      γw    .= γw .+ γl
    end     # g =1:ngs 

    αv,αw = BiorthoScale_vw!(Av,AHw)
    γv[j+1] = αv
    γw[j+1] = αw
  
  else

#   Normalized s.t. <w,v> = 1
    αv,αw = BiorthoScale_vw!(Av,AHw)
    γv[j+1] = αv
    γw[j+1] = αw

  end

end

#----------------------------------------------------------------------
function BiorthoScale_vw!(v::Vector,w::Vector)

    v1      = copy(v)
    w1      = copy(w)

#   Symmetric scaling
    θ       = w'*v                   # <̂w,̂v>
    δ       = sqrt(abs(θ))
    β       = (θ/δ)'

    v      .= v./δ
    w      .= w./β   

    αv      = (v'*v1)/(v'*v)
    αw      = (w'*w1)/(w'*w)

##   ||w|| = 1.0 
#    
#    β       = norm(w)
#    w      .= w./β
#    θ       = w'*v                    # <̂w,̂v>
#    δ       = θ
#    v      .= v./δ
#
#    αv      = δ
#    αw      = β

    return αv,αw

end
#---------------------------------------------------------------------- 






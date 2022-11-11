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

  el  = eltype(u)
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
        Av    .= Av  .- γr[k]*V[:,k]      # ̂v =  A*v_j+1   - β_j*v
        AHw   .= AHw .- γl[k]*W[:,k]      # ̂w = (A')w_j+1  - δ_j*w
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

#  return α,β,δ,u,w

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

  el  = eltype(u)
  zro = el(0.0)

  γv .= zro.*γv
  γw .= zro.*γw

  γr = zro.*copy(γv)
  γl = zro.*copy(γw)

# New Lanczos Vectors
  if j > 0

    for g = 1:ngs
      for k in 1:j
        γr[k]  = W[:,k]'*Av
        γl[k]  = W[:,k]'*AHw
      end
  
      for k in 1:j
        Av    .= Av  .- γr[k]*V[:,k]      # ̂v =  A*v_j+1   - β_j*v
        AHw   .= AHw .- γl[k]*W[:,k]      # ̂w = (A')w_j+1  - δ_j*w
      end 
      γv    .= γv .+ γr
      γw    .= γw .+ γl
    end     # g =1:ngs 

    γv[j+1],γw[j+1] = BiorthoScale_vw!(Av,AHw)
   
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

#    θ       = w'*v                   # <̂w,̂v>
#    β       = norm(w)
#    δ       = θ/(β')

    v      .= v./δ
    w      .= w./β

    αv      = (v'*v1)/(v'*v)
    αw      = (w'*w1)/(w'*w)

    return αv,αw

end
#---------------------------------------------------------------------- 






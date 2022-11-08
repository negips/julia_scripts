# Lanczos Update
function BiOrthoUpd!(Av::Vector,AHw::Vector,V::Matrix,W::Matrix,γv::Vector,γw::Vector,j::Int)

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

  el  = eltype(u)
  zro = el(0.0)

  γv .= zro.*γv
  γw .= zro.*γw

# New Lanczos Vectors
  if j > 0

    for k in 1:j
      γv[k]  = W[:,k]'*Av
      γw[k]  = V[:,k]'*AHw
    end

    for k in 1:j
      Av    .= Av  .- γv[k]*V[:,k]      # ̂v =  A*v_j+1   - β_j*v
      AHw   .= AHw .- γw[k]*W[:,k]      # ̂w = (A')w_j+1  - δ_j*w
    end 
    
    v1      = copy(Av)
    w1      = copy(AHw)

    θ       = AHw'*Av                   # <̂w,̂v>
    δ       = sqrt(abs(θ))
    β       = θ/δ

    Av      .= Av./δ
    AHw     .= AHw./(β')

    γv[j+1] = (Av'*v1)/(Av'*Av)
    γw[j+1] = (AHw'*w1)/(AHw'*AHw)

  else

    θ       = AHw'*Av                   # <̂w,̂v>
    δ       = sqrt(abs(θ))
    β       = θ/δ

#   Normalized s.t. <w,v> = 1
    Av      .= Av./δ
    AHw     .= AHw./(β')
  end

#  return α,β,δ,u,w

end

#----------------------------------------------------------------------



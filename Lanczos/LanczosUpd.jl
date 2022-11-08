# Lanczos Update
function LanczosUpd!(Av::Vector,AHw::Vector,V::Matrix,W::Matrix,γ::Vector,j::Int)

# Yousef Saad - Numerical Methods for large eigenvalue problems (2011) 2nd Edition
# Algorithm 6.6: The non-Hermitian Lanczos Algorithm

# V         - Right Krylov Space
# W         - Left Krylov Space
# Av        - New vector Av     for k>0
# AHw       - New vector (A^H)w for k>0 
# α,β,δ     - Previous/New projections
# j         - Current Krylov size
# jmax      - Maximum Krylov size

  el  = eltype(u)
  zro = el(0.0)
  α   = γ[1]
  δ   = γ[2]
  β   = γ[3]

# New Lanczos Vectors
  if j > 0

    α        = W[:,j]'*Av                # <w,Au>
    if (j > 1) 
      Av     = Av  .- (β )*V[:,j-1]      # ̂v =  A*v_j   - β_j*v_j-1
      AHw    = AHw .- (δ')*W[:,j-1]      # ̂w = (A')w_j  - δ_j*w_j-1
    end  

    Av       = Av  .- (α )*V[:,j]        # ̂v =  A*v_j   - β_j*v_j-1 -  α_j*v_j
    AHw      = AHw .- (α')*W[:,j]        # ̂w = (A')w_j  - δ_j*w_j-1  - (α_j')*w_j

    θ        = AHw'*Av                   # <̂w,̂v>

    δ        = sqrt(abs(θ))
    β        = θ/δ

    Av       = Av./δ
    AHw      = AHw./(β')

  else

    θ        = AHw'*Av                   # <̂w,̂v>

    δ        = sqrt(abs(θ))
    β        = θ/δ

#   Normalized s.t. <w,v> = 1
    Av       = Av./δ
    AHw      = AHw./(β')
    α        = zro
    β        = zro
    δ        = zro
  end

  γ[1] = α
  γ[2] = δ
  γ[3] = β

#  return α,β,δ,u,w

end

#----------------------------------------------------------------------



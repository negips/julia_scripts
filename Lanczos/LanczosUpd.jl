# Lanczos Update
function LanczosUpd!(Q::Matrix,P::Matrix,u::Vector,w::Vector,γ::Vector,k::Int)

# Q         - Right Krylov Space
# P         - Left Krylov Space
# α,β,δ     - Previous/New projections
# k         - Current Krylov size
# kmax      - Maximum Krylov size
# Aq        - New vector u = Aq     for k>0
# ATp       - New vector w = (A^T)p for k>0 

  el  = eltype(u)
  zro = el(0.0)
  α   = γ[1]
  δ   = γ[2]
  β   = γ[3]

# New Lanczos Vectors
  if k > 0

    if (k > 1) 
      u      = u .- δ*Q[:,k-1]
      w      = w .- β*P[:,k-1]
    end  

    α        = P[:,k]'*u
    u        = u .- α*Q[:,k]
    w        = w .- α*P[:,k]

    β        = norm(u)
    δ        = u'*w/β

    u        = u./β
    w        = w./δ

  else
    α        = zro
    β        = norm(u)
    δ        = (u'*w)/β
    
    u        = u./β           # Qs normalized to 1
    w        = w./δ           # Ps normalized s.t. <q,p> = 1

  end
  γ[1] = α
  γ[2] = δ
  γ[3] = β

#  return α,β,δ,u,w

end

#----------------------------------------------------------------------



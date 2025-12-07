# 4th Order Runge-Kutta Steps
# OP is a function acting on v
function OP_RK4!(OP,v,dt)

  localprec = eltype(v[1])
  two = localprec(2)
  six = localprec(6)

  v1 = v .+ dt/two*OP(v)
  v2 = v .+ dt/two*OP(v1)
  v3 = v .+ dt*OP(v2)
  v .= v .+ dt/six*(OP(v) .+ two*OP(v1) .+ two*OP(v2) .+ OP(v3))

  return nothing 
end
#---------------------------------------------------------------------- 
function OP_RK4!(OP,v::AbstractVector{T},dt,wk::AbstractMatrix{T}) where {T <: Number}

  two       = T(2)
  six       = T(6)
  dv1       = view(wk,:,1)
  dv2       = view(wk,:,2)
  dv3       = view(wk,:,3)
  dv4       = view(wk,:,4)
  vn        = view(wk,:,5)

  dv1       = OP(v)
  vn       .= v .+ dt/two*dv1
  dv2       = OP(vn)
  vn       .= v .+ dt/two*dv2
  dv3       = OP(vn)
  vn       .= v .+ dt*dv3
  dv4       = OP(vn)
  v        .= v .+ dt/six*(dv1 .+ two*dv2 .+ two*dv3 .+ dv4)

  return nothing 
end
#---------------------------------------------------------------------- 
function OP_RK4(OP,v::T,dt::Float64) where T <: Union{ComplexF64,Float64}

  two = T(2)
  six = T(6)

  v1 = v + dt/two*OP(v)
  v2 = v + dt/two*OP(v1)
  v3 = v + dt*OP(v2)
  v  = v + dt/six*(OP(v) + two*OP(v1) + two*OP(v2) + OP(v3))

  return v    
end
#----------------------------------------------------------------------
function OP2_RK4!(OP,v::AbstractVector{T},θ::AbstractVector{T},dt::Float64) where T <: Union{ComplexF64,Float64}

  two       = T(2)
  six       = T(6)

  dv,dθ     = OP(v,θ)
  v1        = v + dt/two*dv
  θ1        = θ + dt/two*dθ

  dv1,dθ1   = OP(v1,θ1)
  v2        = v + dt/two*dv1
  θ2        = θ + dt/two*dθ1

  dv2,dθ2   = OP(v2,θ2)
  v3        = v + dt*dv2
  θ3        = θ + dt*dθ2

  dv3,dθ3   = OP(v3,θ3)
  v        .= v .+ dt/six*(dv .+ two*dv1 .+ two*dv2 .+ dv3)
  θ        .= θ .+ dt/six*(dθ .+ two*dθ1 .+ two*dθ2 .+ dθ3)

  return nothing 
end
#----------------------------------------------------------------------
function OP2_RK4!(OP,v::AbstractVector{T},θ::AbstractVector{T},dt::Float64,vw::AbstractMatrix{T},θw::AbstractMatrix{T}) where T <: Union{ComplexF64,Float64}

  two       = T(2)
  six       = T(6)

  dv1 = view(vw,:,1)
  dθ1 = view(θw,:,1)

  dv2 = view(vw,:,2)
  dθ2 = view(θw,:,2)

  dv3 = view(vw,:,3)
  dθ3 = view(θw,:,3)

  dv4 = view(vw,:,4)
  dθ4 = view(θw,:,4)

  vn  = view(vw,:,5)
  θn  = view(θw,:,5)
 
  dv1,dθ1   = OP(v,θ)
  vn        = v + dt/two*dv1
  θn        = θ + dt/two*dθ1

  dv2,dθ2   = OP(vn,θn)
  vn        = v + dt/two*dv2
  θn        = θ + dt/two*dθ2

  dv3,dθ3   = OP(vn,θn)
  vn        = v + dt*dv3
  θn        = θ + dt*dθ3

  dv4,dθ4   = OP(vn,θn)
  v        .= v .+ dt/six*(dv1 .+ two*dv2 .+ two*dv3 .+ dv4)
  θ        .= θ .+ dt/six*(dθ1 .+ two*dθ2 .+ two*dθ3 .+ dθ4)

  return nothing 
end
#----------------------------------------------------------------------





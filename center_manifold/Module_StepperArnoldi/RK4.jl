# 4th Order Runge-Kutta Steps
function RK4!(v::AbstractVector{T},M::AbstractMatrix{T},dt::S) where {T,S<:Number}

  # localprec = eltype(v[1])
  two = T(2)
  six = T(6)

  v1 = v .+ dt/two*M*v
  v2 = v .+ dt/two*M*v1
  v3 = v .+ dt*M*v2
  v .= v .+ dt/six*(M*(v .+ two*v1 .+ two*v2 .+ v3))

  return nothing    
end  
#----------------------------------------------------------------------
function RK4!(v::AbstractVector{T},M::AbstractMatrix{T},v1::AbstractVector{T},v2::AbstractVector{T},v3::AbstractVector{T},dt::S) where {T,S<:Number}

  #localprec = eltype(v[1])
  two = T(2)
  six = T(6)

  v1 .= v .+ dt/two*M*v
  v2 .= v .+ dt/two*M*v1
  v3 .= v .+ dt*M*v2
  v  .= v .+ dt/six*(M*(v .+ two*v1 .+ two*v2 .+ v3))

  return nothing 
end  
#---------------------------------------------------------------------- 
function BiRK4!(v::AbstractVector{T},M::AbstractMatrix{T},Bi::AbstractVector{S},v1::AbstractVector{T},v2::AbstractVector{T},v3::AbstractVector{T},dt::P) where {T,S,P<:Number}

  #localprec = eltype(v[1])
  two = T(2)
  six = T(6)

  v1 .= v .+ dt/two*Bi.*(M*v)
  v2 .= v .+ dt/two*Bi.*(M*v1)
  v3 .= v .+ dt*Bi.*(M*v2)
  v  .= v .+ dt/six*Bi.*(M*(v .+ two*v1 .+ two*v2 .+ v3))

  return nothing 
end  
#---------------------------------------------------------------------- 
function RestrictedBRK4!(v::AbstractVector{T},M::AbstractMatrix{T},B::AbstractVector{S},Vr::AbstractMatrix{T},Wr::AbstractMatrix{T},lbc::Bool,rbc::Bool,v1::AbstractVector{T},v2::AbstractVector{T},v3::AbstractVector{T},dt::P) where {T,S,P<:Number}

  #localprec = eltype(v[1])
  two = T(2)
  six = T(6)
  ngs = 2

  Bi  = 1.0./B

  v1 .= v .+ dt/two*(Bi.*(M*v))
  ObliqueSubspaceRemoval!(v1,Vr,Wr,B,ngs) 
  StpArn_SetBC!(v1,lbc,rbc)

  v2 .= v .+ dt/two*Bi.*(M*v1)
  ObliqueSubspaceRemoval!(v2,Vr,Wr,B,ngs) 
  StpArn_SetBC!(v2,lbc,rbc)
 
  v3 .= v .+ dt*Bi.*(M*v2)
  ObliqueSubspaceRemoval!(v3,Vr,Wr,B,ngs) 
  StpArn_SetBC!(v3,lbc,rbc)
 
  v  .= v .+ dt/six*Bi.*(M*(v .+ two*v1 .+ two*v2 .+ v3))
  ObliqueSubspaceRemoval!(v,Vr,Wr,B,ngs) 
  StpArn_SetBC!(v,lbc,rbc)

  return nothing 
end  
#---------------------------------------------------------------------- 










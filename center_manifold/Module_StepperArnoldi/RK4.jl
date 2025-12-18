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
function RK4!(v::AbstractVector{T},OP,dt::S) where {T,S<:Number}

  # localprec = eltype(v[1])
  two = T(2)
  six = T(6)

  v1 = v .+ dt/two*OP(v)
  v2 = v .+ dt/two*OP(v1)
  v3 = v .+ dt*OP(v2)
  v .= v .+ dt/six*(OP(v) .+ two*OP(v1) .+ two*OP(v2) .+ OP(v3))

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
function RK4!(v::AbstractVector{T},OP,v1::AbstractVector{T},v2::AbstractVector{T},v3::AbstractVector{T},dt::S) where {T,S<:Number}

  #localprec = eltype(v[1])
  two = T(2)
  six = T(6)

  v1 .= v .+ dt/two*OP(v)
  v2 .= v .+ dt/two*OP(v1)
  v3 .= v .+ dt*OP(v2)
  v  .= v .+ dt/six*(OP(v) .+ two*OP(v1) .+ two*OP(v2) .+ OP(v3))

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
function BiRK4!(v::AbstractVector{T},OP,Bi::AbstractVector{S},v1::AbstractVector{T},v2::AbstractVector{T},v3::AbstractVector{T},dt::P) where {T,S,P<:Number}

  #localprec = eltype(v[1])
  two = T(2)
  six = T(6)

  v1 .= v .+ dt/two*Bi.*OP(v)
  v2 .= v .+ dt/two*Bi.*OP(v1)
  v3 .= v .+ dt*Bi.*OP(v2)
  v  .= v .+ dt/six*Bi.*(OP(v) .+ two*OP(v1) .+ two*OP(v2) .+ OP(v3))

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
function RestrictedBRK4!(v::AbstractVector{T},OP,B::AbstractVector{S},Vr::AbstractMatrix{T},Wr::AbstractMatrix{T},lbc::Bool,rbc::Bool,v1::AbstractVector{T},v2::AbstractVector{T},v3::AbstractVector{T},dt::P) where {T,S,P<:Number}

  #localprec = eltype(v[1])
  two = T(2)
  six = T(6)
  ngs = 2

  Bi  = 1.0./B

  v1 .= v .+ dt/two*Bi.*OP(v)
  ObliqueSubspaceRemoval!(v1,Vr,Wr,B,ngs) 
  StpArn_SetBC!(v1,lbc,rbc)

  v2 .= v .+ dt/two*Bi.*OP(v1)
  ObliqueSubspaceRemoval!(v2,Vr,Wr,B,ngs) 
  StpArn_SetBC!(v2,lbc,rbc)
 
  v3 .= v .+ dt*Bi.*OP(v2)
  ObliqueSubspaceRemoval!(v3,Vr,Wr,B,ngs) 
  StpArn_SetBC!(v3,lbc,rbc)
 
  v  .= v .+ dt/six*Bi.*(OP(v) .+ two*OP(v1) .+ two*OP(v2) .+ OP(v3))
  ObliqueSubspaceRemoval!(v,Vr,Wr,B,ngs) 
  StpArn_SetBC!(v,lbc,rbc)

  return nothing 
end  
#---------------------------------------------------------------------- 
function REP_BRK4!(ve::AbstractVector{T1},L::AbstractMatrix{T2},B::AbstractVector{T3},σ::AbstractVector{T1},V::AbstractMatrix{T1},W::AbstractMatrix{T1},restriction::Vector{Bool},f::AbstractVector{T2},ω::T1,lbc::Bool,rbc::Bool,v1::AbstractVector{T1},v2::AbstractVector{T1},v3::AbstractVector{T1},v4::AbstractVector{T1},v5::AbstractVector{T1},dt::P) where {T,S,P<:Number}

  Ne  = length(ve)
  N   = Ne - 1
  two = T(2)
  six = T(6)
  ngs = 2

  Bi  = 1.0./B

  v1  = REPLx(ve,L,B,V,W,σ,restricted,f,ω)
  for i in 1:N
    v1[i] = Bi[i]*v1[i]
    v5[i] = ve[i] + dt/two*v1[i]
  end
  v5[Ne]  = ve[Ne] + dt/two*v1[Ne]
  @views ObliqueSubspaceRemoval!(v5[1:N],V,W,B,ngs) 
  @views StpArn_SetBC!(v5[1:N],lbc,rbc)

  v2  = REPLx(v5,L,B,V,W,σ,restricted,f,ω)
  for i in 1:N
    v2[i] = Bi[i]*v2[i]
    v5[i] = ve[i] + dt/two*v2[i]
  end
  v5[Ne]  = ve[Ne] + dt/two*v2[Ne]
  @views ObliqueSubspaceRemoval!(v5[1:N],V,W,B,ngs) 
  @views StpArn_SetBC!(v5[1:N],lbc,rbc)

  v3  = REPLx(v5,L,B,V,W,σ,restricted,f,ω)
  for i in 1:N
    v3[i] = Bi[i]*v3[i]
    v5[i] = ve[i] + dt*v3[i]
  end
  v5[Ne]  = ve[Ne] + dt*v3[Ne]
  @views ObliqueSubspaceRemoval!(v5[1:N],V,W,B,ngs) 
  @views StpArn_SetBC!(v5[1:N],lbc,rbc)

  v4  = REPLx(v5,L,B,V,W,σ,restricted,f,ω)
  for i in 1:N
    v4[i] = Bi[i]*v4[i]
  end

  for i in 1:Ne
    v[i]  .= v[i] .+ dt/six*(v1 .+ two*v2 .+ two*v3 .+ v4)
  end  
  @views ObliqueSubspaceRemoval!(v[1:N],V,W,B,ngs) 
  @views StpArn_SetBC!(v[1:N],lbc,rbc)

  return nothing 
end  
#---------------------------------------------------------------------- 









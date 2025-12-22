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
function BRK4_2!(ve::AbstractVector{T1},M::AbstractMatrix{T2},Be::AbstractVector{T3},f::AbstractVector{T2},ω::T2,ve1::AbstractVector{T1},ve2::AbstractVector{T1},ve3::AbstractVector{T1},dt::T4) where {T1,T2,T3,T4<:Number}

  #localprec = eltype(v[1])
  two = T1(2)
  six = T1(6)

  N   = size(M,2)
  Ne  = N+1

  v    = view(ve,1:N)
  v1   = view(ve1,1:N)
  v2   = view(ve2,1:N)
  v3   = view(ve3,1:N)

  Bei = 1.0./Be
  B   = view(Be,1:N)
  Bi  = view(Bei,1:N)

  v1       .= v .+ dt/two*Bi.*(M*v .+ B.*f*ve[Ne])
  ve1[Ne]   = ve[Ne] + dt/two*ω*ve[Ne]

  v2       .= v .+ dt/two*Bi.*(M*v1 .+ B.*f*ve1[Ne])
  ve2[Ne]   = ve[Ne] + dt/two*ω*ve1[Ne]

  v3       .= v .+ dt*Bi.*(M*v2 .+ B.*f*ve2[Ne])
  ve3[Ne]   = ve[Ne] + dt*ω*ve2[Ne]

  v        .= v .+ dt/six*Bi.*(M*(v .+ two*v1 .+ two*v2 .+ v3) .+ B.*f*(ve[Ne] .+ two*ve1[Ne] .+ two*ve2[Ne] .+ ve3[Ne]))
  ve[Ne]    = ve[Ne] + dt/six*ω*(ve[Ne] + two*ve1[Ne] + two*ve2[Ne] + ve3[Ne])

  return nothing 
end  
#---------------------------------------------------------------------- 
function BRK4_3!(ve::AbstractVector{T1},M::AbstractMatrix{T2},Be::AbstractVector{T3},f::AbstractVector{T2},ω::T2,ve1::AbstractVector{T1},ve2::AbstractVector{T1},ve3::AbstractVector{T1},ve4::AbstractVector{T1},dt::T4) where {T1,T2,T3,T4<:Number}

  #localprec = eltype(v[1])
  two = T1(2)
  six = T1(6)

  N   = size(M,2)
  Ne  = N+1

  v    = view(ve,1:N)
  v1   = view(ve1,1:N)
  v2   = view(ve2,1:N)
  v3   = view(ve3,1:N)
  v4   = view(ve4,1:N)

  Bei = 1.0./Be
  B   = view(Be,1:N)
  Bi  = view(Bei,1:N)

  #mul!(v4,M,v)
  PLx!(v4,v,M,B,zeros(T1,N,0),zeros(T1,N,0),zeros(T1,0))
  v1       .= v .+ dt/two*Bi.*(v4 .+ B.*f*ve[Ne])
  ve1[Ne]   = ve[Ne] + dt/two*ω*ve[Ne]

  #mul!(v4,M,v1)
  PLx!(v4,v1,M,B,zeros(T1,N,0),zeros(T1,N,0),zeros(T1,0))
  v2       .= v .+ dt/two*Bi.*(v4 .+ B.*f*ve1[Ne])
  ve2[Ne]   = ve[Ne] + dt/two*ω*ve1[Ne]

  # mul!(v4,M,v2)
  PLx!(v4,v2,M,B,zeros(T1,N,0),zeros(T1,N,0),zeros(T1,0))
  v3       .= v .+ dt*Bi.*(v4 .+ B.*f*ve2[Ne])
  ve3[Ne]   = ve[Ne] + dt*ω*ve2[Ne]


  ve3 .= ve .+ two*ve1 .+ two*ve2 .+ ve3
  #mul!(v4,M,v3)
  PLx!(v4,v3,M,B,zeros(T1,N,0),zeros(T1,N,0),zeros(T1,0))

  v        .= v .+ dt/six*Bi.*(v4 .+ B.*f*(ve3[Ne]))
  ve[Ne]    = ve[Ne] + dt/six*ω*(ve3[Ne])

  return nothing 
end  
#---------------------------------------------------------------------- 
function RE_BRK4!(ve::AbstractVector{T1},L::AbstractMatrix{T2},Be::AbstractVector{T3},V::AbstractMatrix{T1},W::AbstractMatrix{T1},restricted::Vector{Bool},f::AbstractVector{T2},ω::T1,lbc::Bool,rbc::Bool,ve1::AbstractVector{T1},ve2::AbstractVector{T1},ve3::AbstractVector{T1},ve4::AbstractVector{T1},dt::T4) where {T1,T2,T3,T4<:Number}

  two = T1(2)
  six = T1(6)
  ngs = 2

  N   = size(L,2)
  Ne  = N+1

  v    = view(ve,1:N)
  v1   = view(ve1,1:N)
  v2   = view(ve2,1:N)
  v3   = view(ve3,1:N)
  v4   = view(ve4,1:N)

  Bei = 1.0./Be
  B   = view(Be,1:N)
  Bi  = view(Bei,1:N)

  # Perturbed Lv
  PLx!(v4,v,L,B,V,W,σ)
  v1       .= v .+ dt/two*Bi.*(v4 .+ B.*f*ve[Ne])
  ve1[Ne]   = ve[Ne] + dt/two*ω*ve[Ne]
  ObliqueSubspaceRemoval2!(v1,V,W,B,restricted,ngs) 

  # Perturbed Lv
  PLx!(v4,v1,L,B,V,W,σ)
  v2       .= v .+ dt/two*Bi.*(v4 .+ B.*f*ve1[Ne])
  ve2[Ne]   = ve[Ne] + dt/two*ω*ve1[Ne]
  ObliqueSubspaceRemoval2!(v2,V,W,B,restricted,ngs) 

  # Perturbed Lv
  PLx!(v4,v2,L,B,V,W,σ)
  v3       .= v .+ dt*Bi.*(v4 .+ B.*f*ve2[Ne])
  ve3[Ne]   = ve[Ne] + dt*ω*ve2[Ne]
  ObliqueSubspaceRemoval2!(v3,V,W,B,restricted,ngs) 

  # Perturbed Lv
  ve3 .= ve .+ two*ve1 .+ two*ve2 .+ ve3
  PLx!(v4,v3,L,B,V,W,σ)
  v        .= v .+ dt/six*Bi.*(v4 .+ B.*f*ve3[Ne])
  ve[Ne]    = ve[Ne] + dt/six*ω*(ve3[Ne])
  ObliqueSubspaceRemoval2!(v,V,W,B,restricted,ngs) 

  return nothing 

end  
#---------------------------------------------------------------------- 
function E_BRK4!(ve::AbstractVector{T1},L::AbstractMatrix{T2},Be::AbstractVector{T3},f::AbstractVector{T2},ω::T1,lbc::Bool,rbc::Bool,ve1::AbstractVector{T1},ve2::AbstractVector{T1},ve3::AbstractVector{T1},ve4::AbstractVector{T1},dt::T4) where {T1,T2,T3,T4<:Number}

  two = T1(2)
  six = T1(6)
  ngs = 2

  N   = size(L,2)
  Ne  = N+1

  Bei = 1.0./Be
  B   = view(Be,1:N)
  Bi  = view(Bei,1:N)

  # Extended Lv
  ELx!(ve4,ve,L,B,f,ω)
  for j in 1:N
    ve1[j]  = ve[j] + dt/two*Bei[j]*ve4[j]
  end
  ve1[Ne]  = ve[Ne] + dt/two*ve4[Ne]

  # Extended Lv
  ELx!(ve4,ve1,L,B,f,ω)
  for j in 1:N
    ve2[j]  = ve[j] + dt/two*Bei[j]*ve4[j]
  end
  ve2[Ne]  = ve[Ne] + dt/two*ve4[Ne]

  # Extended Lv
  ELx!(ve4,ve2,L,B,f,ω)
  for j in 1:N
    ve3[j]  = ve[j] + dt*Bei[j]*ve4[j]
  end
  ve3[Ne]  = ve[Ne] + dt*ve4[Ne]

  # Extended Lv
  ve3 .= ve .+ two*ve1 .+ two*ve2 .+ ve3
  ELx!(ve4,ve3,L,B,f,ω)
  for j in 1:N
    ve[j]  = ve[j] + dt/six*Bei[j]*ve4[j]
  end
  ve[Ne]  = ve[Ne] + dt/six*ve4[Ne]

  return nothing 
end  
#---------------------------------------------------------------------- 
function RE_BRK4!(ve::AbstractVector{T1},L::AbstractMatrix{T2},Be::AbstractVector{T3},V::AbstractMatrix{T1},W::AbstractMatrix{T1},restricted::Vector{Bool},f::AbstractVector{T2},ω::T1,lbc::Bool,rbc::Bool,ve1::AbstractVector{T1},ve2::AbstractVector{T1},ve3::AbstractVector{T1},ve4::AbstractVector{T1},ve5::AbstractVector{T1},dt::T4) where {T1,T2,T3,T4<:Number}

  two = T1(2)
  six = T1(6)
  ngs = 2

  N   = size(L,2)
  Ne  = N+1

  Bei = 1.0./Be
  B   = view(Be,1:N)
  Bi  = view(Bei,1:N)

  v    = view(ve ,1:N)
  v1   = view(ve1,1:N)
  v2   = view(ve2,1:N)
  v3   = view(ve3,1:N)
  v4   = view(ve4,1:N)

  # ve5 used as work array

  # Restricted Extended Lv
  RELx!(ve4,ve,L,B,V,W,restricted,f,ω,ve5)
  for j in 1:N
    ve1[j]  = ve[j] + dt/two*Bei[j]*ve4[j]
  end
  ve1[Ne]  = ve[Ne] + dt/two*ve4[Ne]
  ObliqueSubspaceRemoval2!(v1,V,W,B,restricted,ngs)

  # Restricted Extended Lv
  RELx!(ve4,ve1,L,B,V,W,restricted,f,ω,ve5)
  for j in 1:N
    ve2[j]  = ve[j] + dt/two*Bei[j]*ve4[j]
  end
  ve2[Ne]  = ve[Ne] + dt/two*ve4[Ne]
  ObliqueSubspaceRemoval2!(v2,V,W,B,restricted,ngs)

  # Restricted Extended Lv
  RELx!(ve4,ve2,L,B,V,W,restricted,f,ω,ve5)
  for j in 1:N
    ve3[j]  = ve[j] + dt*Bei[j]*ve4[j]
  end
  ve3[Ne]  = ve[Ne] + dt*ve4[Ne]
  ObliqueSubspaceRemoval2!(v3,V,W,B,restricted,ngs)

  # Restricted Extended Lv
  ve3 .= ve .+ two*ve1 .+ two*ve2 .+ ve3
  RELx!(ve4,ve3,L,B,V,W,restricted,f,ω,ve5)
  for j in 1:N
    ve[j]  = ve[j] + dt/six*Bei[j]*ve4[j]
  end
  ve[Ne]  = ve[Ne] + dt/six*ve4[Ne]
  ObliqueSubspaceRemoval2!(v,V,W,B,restricted,ngs)
 

  return nothing 
end  
#---------------------------------------------------------------------- 
function REP_BRK4!(ve::AbstractVector{T1},L::AbstractMatrix{T2},Be::AbstractVector{T3},σ::AbstractVector{T1},V::AbstractMatrix{T1},W::AbstractMatrix{T1},restricted::Vector{Bool},f::AbstractVector{T2},ω::T1,lbc::Bool,rbc::Bool,ve1::AbstractVector{T1},ve2::AbstractVector{T1},ve3::AbstractVector{T1},ve4::AbstractVector{T1},dt::T4) where {T1,T2,T3,T4<:Number}

  two = T1(2)
  six = T1(6)
  ngs = 2

  N   = size(L,2)
  Ne  = N+1

  v    = view(ve,1:N)
  v1   = view(ve1,1:N)
  v2   = view(ve2,1:N)
  v3   = view(ve3,1:N)
  v4   = view(ve4,1:N)

  Bei = 1.0./Be
  B   = view(Be,1:N)
  Bi  = view(Bei,1:N)

  # Perturbed Lv
  PLx!(v4,v,L,B,V,W,σ)
  v1       .= v .+ dt/two*Bi.*(v4 .+ B.*f*ve[Ne])
  ve1[Ne]   = ve[Ne] + dt/two*ω*ve[Ne]
  ObliqueSubspaceRemoval2!(v1,V,W,B,restricted,ngs) 

  # Perturbed Lv
  PLx!(v4,v1,L,B,V,W,σ)
  v2       .= v .+ dt/two*Bi.*(v4 .+ B.*f*ve1[Ne])
  ve2[Ne]   = ve[Ne] + dt/two*ω*ve1[Ne]
  ObliqueSubspaceRemoval2!(v2,V,W,B,restricted,ngs) 

  # Perturbed Lv
  PLx!(v4,v2,L,B,V,W,σ)
  v3       .= v .+ dt*Bi.*(v4 .+ B.*f*ve2[Ne])
  ve3[Ne]   = ve[Ne] + dt*ω*ve2[Ne]
  ObliqueSubspaceRemoval2!(v3,V,W,B,restricted,ngs) 

  # Perturbed Lv
  ve3 .= ve .+ two*ve1 .+ two*ve2 .+ ve3
  PLx!(v4,v3,L,B,V,W,σ)
  v        .= v .+ dt/six*Bi.*(v4 .+ B.*f*ve3[Ne])
  ve[Ne]    = ve[Ne] + dt/six*ω*(ve3[Ne])
  ObliqueSubspaceRemoval2!(v,V,W,B,restricted,ngs) 

  return nothing 

end  
#---------------------------------------------------------------------- 









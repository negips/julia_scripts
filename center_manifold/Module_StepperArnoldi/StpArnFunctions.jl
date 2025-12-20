# Function Definitions
#---------------------------------------------------------------------- 
function copy(Ai::ArnoldiInput)

  Ao = ArnoldiInput(Ai.ifarnoldi,
                    Ai.ifverbose,
                    Ai.ifeigshift,
                    Ai.vlen,
                    Ai.nev,
                    Ai.ekryl,
                    Ai.lkryl,
                    Ai.eigshift,
                    Ai.ngs,
                    Ai.bsize,
                    Ai.outer_iterations,
                    Ai.tol)

  return Ao
end
#----------------------------------------------------------------------
function copy(Si::StepperInput)

  So = StepperInput(Si.ifadjoint,
                    Si.ifoptimal,
                    Si.ifverbose,
                    Si.verbosestep,
                    Si.nsteps,
                    Si.timestep)

  return So
end  
#---------------------------------------------------------------------- 
function ObliqueSubspaceRemoval!(v::AbstractVector{T},V::AbstractMatrix{T},W::AbstractMatrix{T},B::AbstractVector{S},ngs::Int) where {T,S<:Number}

  for i in 1:ngs
    α   = W'*(B.*v)
    v  .= v .- V*α
  end  

  return nothing
end  
#---------------------------------------------------------------------- 
function ObliqueSubspaceRemoval2!(v::AbstractVector{T},V::AbstractMatrix{T},W::AbstractMatrix{T},B::AbstractVector{S},ifremove::Vector{Bool},ngs::Int) where {T,S<:Number}

  nv = size(V,2)

  for j in 1:ngs
    β    = zeros(T,nv)
    for i in 1:nv
      if ifremove[i]
        β[i]  = W[:,i]'*(B.*v)
      end
    end   

    for i in 1:nv
      if ifremove[i]
        v  .= v .- V[:,i]*β[i]
      end
    end
  end

  return nothing
end  
#---------------------------------------------------------------------- 
function ObliqueSubspaceRemoval3!(v::AbstractVector{T},V::AbstractMatrix{T},W::AbstractMatrix{T},B::AbstractVector{S},ifremove::Vector{Bool},ngs::Int) where {T,S<:Number}

  nv = size(V,2)

  # Get components of subspace V in v
  α    = zeros(T,nv)
  for j in 1:ngs
    β    = zeros(T,nv)
    for i in 1:nv
      if ifremove[i]
        β[i]  = W[:,i]'*(B.*v)
        α[i]  = α[i] + β[i]
      end
    end   

    for i in 1:nv
      if ifremove[i]
        v  .= v .- V[:,i]*β[i]
      end
    end
  end

  return α
end  
#---------------------------------------------------------------------- 
function PLx(x::AbstractVector{T1},L::AbstractMatrix{T1},B::AbstractVector{T2},V::AbstractMatrix{T1},W::AbstractMatrix{T1},σ::AbstractVector{T3}) where {T1,T2,T3<:Number}

  ngs  = 2
  N    = size(L,1)
  nv   = size(V,2)
  σtol = 1.0e-12

  # using Lx as a work array for now
  Lx   = copy(x)

  # Get components of the Perturbed Modes in x
  α    = zeros(T1,nv)
  for j in 1:ngs
    β    = zeros(T1,nv)
    for i in 1:nv
      if abs(σ[i]) > σtol
        β[i]  = W[:,i]'*(B.*Lx)
        α[i]  = α[i] + β[i]
      end
    end   

    # Since we only need α.
    if j<ngs
      for i in 1:nv
        if abs(σ[i]) > σtol
          Lx  .= Lx .- V[:,i]*β[i]
        end
      end
    end  
  end
  # Lx = L*x
  mul!(Lx,L,x)

  # Add Mode Perturbations
  for i in 1:nv
    if abs(σ[i]) > σtol
      Lx  .= Lx .+ σ[i]*B.*(V[:,i]*α[i])
    end
  end   

  return Lx
end
#---------------------------------------------------------------------- 
function PLx!(Lx::AbstractVector{T1},x::AbstractVector{T1},L::AbstractMatrix{T1},B::AbstractVector{T2},V::AbstractMatrix{T1},W::AbstractMatrix{T1},σ::AbstractVector{T3}) where {T1,T2,T3<:Number}

  ngs  = 2
  N    = size(L,1)
  nv   = size(V,2)
  σtol = 1.0e-12

  # using Lx as a work array for now
  copyto!(Lx,1,x,1,N)

  # Get components of the Perturbed Modes in x
  α    = zeros(T1,nv)
  for j in 1:ngs
    β    = zeros(T1,nv)
    for i in 1:nv
      if abs(σ[i]) > σtol
        β[i]  = W[:,i]'*(B.*Lx)
        α[i]  = α[i] + β[i]
      end
    end   

    # Since we only need α.
    if j<ngs
      for i in 1:nv
        if abs(σ[i]) > σtol
          Lx  .= Lx .- V[:,i]*β[i]
        end
      end
    end  
  end

  # Lx = L*x
  mul!(Lx,L,x)

  # Add Mode Perturbations
  for i in 1:nv
    if abs(σ[i]) > σtol
      Lx  .= Lx .+ σ[i]*B.*(V[:,i]*α[i])
    end
  end   

  return nothing
end
#---------------------------------------------------------------------- 
function RPLx!(Lx::AbstractVector{T1},x::AbstractVector{T1},L::AbstractMatrix{T1},B::AbstractVector{T2},V::AbstractMatrix{T1},W::AbstractMatrix{T1},σ::AbstractVector{T3},restricted::Vector{Bool},x1::AbstractVector{T1}) where {T1,T2,T3<:Number}

  ngs  = 2
  N    = size(L,1)
  nv   = size(V,2)
  σtol = 1.0e-12

  # using x1 as a work array for now
  copyto!(x1,1,x,1,N)

  # Get components of the Perturbed Modes in x
  α    = zeros(T1,nv)
  for j in 1:ngs
    β    = zeros(T1,nv)
    for i in 1:nv
      if abs(σ[i]) > σtol && ~restricted[i]
        β[i]  = W[:,i]'*(B.*x1)
        α[i]  = α[i] + β[i]
      end
    end   

    # Since we only need α.
    if j<ngs
      for i in 1:nv
        if abs(σ[i]) > σtol
          x1  .= x1 .- V[:,i]*β[i]
        end
      end
    end  
  end

  copyto!(x1,1,x,1,N)
  ObliqueSubspaceRemoval2!(x1,V,W,B,restriction,ngs)

  # Lx = L*x
  mul!(Lx,L,x1)

  # Add Mode Perturbations
  for i in 1:nv
    if abs(σ[i]) > σtol
      Lx  .= Lx .+ σ[i]*B.*(V[:,i]*α[i])
    end
  end   

  return nothing
end
#---------------------------------------------------------------------- 
function EPLx(xe::AbstractVector{T1},L::AbstractMatrix{T1},B::AbstractVector{T2},V::AbstractMatrix{T1},W::AbstractMatrix{T1},σ::AbstractVector{T3},f::AbstractVector{T1},ω::T3) where {T1,T2,T3<:Number}

  ngs = 2
  N   = size(L,2)
  Ne  = N+1
  Lx  = zeros(T1,Ne)
  np  = length(σ)
  σtol = 1.0e-12
  
  @views PLx!(Lx[1:N],xe[1:N],L,B,V,W,σ)

  # Add forcing extension
  for j in 1:N
    Lx[j] = Lx[j] + B[j]*f[j]*xe[Ne]
  end
  Lx[Ne]  = ω*xe[Ne]

  return Lx
end
#---------------------------------------------------------------------- 
function EPLx!(Lxe::AbstractVector{T1},xe::AbstractVector{T1},L::AbstractMatrix{T1},B::AbstractVector{T2},V::AbstractMatrix{T1},W::AbstractMatrix{T1},σ::AbstractVector{T3},f::AbstractVector{T1},ω::T3) where {T1,T2,T3<:Number}

  ngs  = 2
  N    = size(L,2)
  Ne   = N+1
  Lxe  = zeros(T1,Ne)
  nv   = size(V,2)
  σtol = 1.0e-12
  
  @views PLx!(Lxe[1:N],xe[1:N],L,B,V,W,σ)

  # Add forcing extension
  for j in 1:N
    Lxe[j] = Lxe[j] + B[j]*f[j]*xe[Ne]
  end
  Lxe[Ne]  = ω*xe[Ne]

  return nothing
end
#---------------------------------------------------------------------- 

function REPLx(xe::AbstractVector{T1},L::AbstractMatrix{T1},B::AbstractVector{T2},V::AbstractMatrix{T1},W::AbstractMatrix{T1},σ::AbstractVector{T3},restricted::AbstractVector{Bool},f::AbstractVector{T1},ω::T3) where {T1,T2,T3<:Number}

  ngs = 2
  N   = size(L,2)
  Ne  = N+1
  Lx  = zeros(T1,Ne)
  σtol = 1.0e-12
  
  @views RPLx!(Lx[1:N],xe[1:N],L,B,V,W,σ,restricted)

  # Add forcing extension
  # ftmp = copy(f)
  # ObliqueSubspaceRemoval2!(ftmp,V,W,B,restricted,ngs)
  for j in 1:N
    Lx[j] = Lx[j] + B[j]*f[j]*xe[Ne]
  end
  Lx[Ne]  = ω*xe[Ne]

  return Lx
end
#---------------------------------------------------------------------- 













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
    for i in 1:nv
      if (ifremove[i])
        β   = W[:,i]'*(B.*v)
        v  .= v .- V[:,i]*β
      end
    end
  end

  return nothing
end  
#---------------------------------------------------------------------- 
function ObliqueSubspaceRemoval3!(v::AbstractVector{T},V::AbstractMatrix{T},W::AbstractMatrix{T},B::AbstractVector{S},ifremove::Vector{Bool},ngs::Int) where {T,S<:Number}

  nv = size(V,2)

  α  = zeros(T,nv)

  for j in 1:ngs
    for i in 1:nv
      if (ifremove[i])
        β     = W[:,i]'*(B.*v)
        v    .= v .- V[:,i]*β
        α[i]  = α[i] + β
      end
    end
  end

  return α
end  
#---------------------------------------------------------------------- 
function PLx(x::AbstractVector{T1},L::AbstractMatrix{T1},B::AbstractVector{T2},V::AbstractMatrix{T1},W::AbstractMatrix{T1},σ::AbstractVector{T3}) where {T1,T2,T3<:Number}

  ngs = 2
  N   = size(L,2)
  Lx  = zeros(T1,N)
  np  = length(σ)
  σtol = 1.0e-12
  xtmp = Vector{T1}(undef,N)

  # L*x
  Lx = L*x

  # Evaluate Mode perturbation terms
  for i in 1:np
    if abs(σ[i]) > σtol
      copyto!(xtmp,1,x,1,N)
      α    = T1(0)
      for j in 1:ngs
        β     = W[:,i]'*(B.*xtmp)
        α     = α + β 
        xtmp .= xtmp .- β*V[:,i]
      end
      for j in 1:N
        Lx[j] = Lx[j] + α*σ[i]*B[j]*V[j,i]
      end
    end
  end  

  return Lx
end
#---------------------------------------------------------------------- 
function PLx!(Lx::AbstractVector{T1},x::AbstractVector{T1},L::AbstractMatrix{T1},B::AbstractVector{T2},V::AbstractMatrix{T1},W::AbstractMatrix{T1},σ::AbstractVector{T3}) where {T1,T2,T3<:Number}

  ngs  = 2
  N    = size(L,2)
  nv   = size(V,2)
  σtol = 1.0e-12

  # Lx = L*x
  mul!(Lx,L,x)

  # Evaluate Mode perturbation terms
  xtmp = Vector{T1}(undef,N)
  for i in 1:nv
    if abs(σ[i]) > σtol
      copyto!(xtmp,1,x,1,N)
      α    = T1(0)
      for j in 1:ngs
        β     = W[:,i]'*(B.*xtmp)
        α     = α + β 
        xtmp .= xtmp .- β*V[:,i]
      end
      for j in 1:N
        Lx[j] = Lx[j] + α*σ[i]*B[j]*V[j,i]
      end
    end
  end  

  return nothing
end
#---------------------------------------------------------------------- 
function RPLx!(Lx::AbstractVector{T1},x::AbstractVector{T1},L::AbstractMatrix{T1},B::AbstractVector{T2},V::AbstractMatrix{T1},W::AbstractMatrix{T1},σ::AbstractVector{T3},restricted::Vector{Bool}) where {T1,T2,T3<:Number}

  ngs  = 2
  N    = size(L,2)
  np   = length(σ)
  σtol = 1.0e-12

  # Lx = L*x
  xtmp = Vector{T1}(undef,N)
  copyto!(xtmp,1,x,1,N)

  ObliqueSubspaceRemoval2!(xtmp,V,W,B,restricted,ngs)
  mul!(Lx,L,xtmp)

  # Evaluate Mode perturbation terms
  for i in 1:np
    if (~restricted[i])
      if abs(σ[i]) > σtol
        copyto!(xtmp,1,x,1,N)
        α    = T1(0)
        for j in 1:ngs
          β     = W[:,i]'*(B.*xtmp)
          α     = α + β 
          xtmp .= xtmp .- β*V[:,i]
        end
        for j in 1:N
          Lx[j] = Lx[j] + α*σ[i]*B[j]*V[j,i]
        end
      end         # abs(σ[i]) > σtol
    end           # ~restricted[i]
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
  
  @views PertLx!(Lx[1:N],xe[1:N],L,B,V,W,σ)

  # Add forcing extension
  for j in 1:N
    Lx[j] = Lx[j] + B[j]*f[j]
  end
  Lx[Ne]  = ω*x[Ne]

  return Lx
end
#---------------------------------------------------------------------- 
function REPLx(xe::AbstractVector{T1},L::AbstractMatrix{T1},B::AbstractVector{T2},V::AbstractMatrix{T1},W::AbstractMatrix{T1},σ::AbstractVector{T3},restricted::AbstractVector{Bool},f::AbstractVector{T1},ω::T3) where {T1,T2,T3<:Number}

  ngs = 2
  N   = size(L,2)
  Ne  = N+1
  Lx  = zeros(T1,Ne)
  nv  = size(V,2)
  σtol = 1.0e-12
  
  @views RPLx!(Lx[1:N],xe[1:N],L,B,V,W,σ,restricted)

  # Add forcing extension
  ftmp = copy(f)
  ObliqueSubspaceRemoval2!(ftmp,V,W,B,restricted,ngs)
  for j in 1:N
    Lx[j] = Lx[j] + B[j]*ftmp[j]
  end
  Lx[Ne]  = ω*x[Ne]

  return Lx
end
#---------------------------------------------------------------------- 













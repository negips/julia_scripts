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













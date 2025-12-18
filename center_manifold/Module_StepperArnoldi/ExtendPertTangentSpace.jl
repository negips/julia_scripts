# Extending the tangent space
#---------------------------------------------------------------------- 
function EPTangentSpaces(L::AbstractMatrix{T2},B::AbstractVector{T3},λc::AbstractVector{T1},σ::AbstractVector{T1},V::AbstractMatrix{T2},W::AbstractMatrix{T2},F::AbstractMatrix{T2},λe::AbstractVector{T1},restricted::Bool,AInp::ArnoldiInput,SInp::StepperInput,lbc::Bool,rbc::Bool) where {T1,T2,T3<:Number}

  ne = length(λe)

  EM = Vector{StepperArnoldi.ExtendedMode}(undef,ne)

  for i in 1:ne
    f     = view(F,:,i)
    AInp.ifeigshift = true
    AInp.eigshift   = λe[i]
    if (restricted)
      EM[i] = ExtendTangentSpaceRestricted(L,B,λc,V,W,f,λe[i],AInp,SInp,lbc,rbc)
    else
      EM[i] = EPTangentSpace(L,B,λc,σ,V,W,f,λe[i],AInp,SInp,lbc,rbc)
    end
  end  

  EModes = ExtendedModes(EM)

  return EModes
end
#---------------------------------------------------------------------- 
function EPTangentSpace(L::AbstractMatrix{T2},B::AbstractVector{T3},λc::AbstractVector{T1},σ::AbstractVector{T1},V::AbstractMatrix{T2},W::AbstractMatrix{T2},f::AbstractVector{T2},λe::T1,AInp::ArnoldiInput,SInp::StepperInput,lbc::Bool,rbc::Bool) where {T1,T2,T3<:Number}

  ngs = 2
  N   = size(L,2)
  n   = length(λc)

  Γ   = zeros(T2,n)
  vh  = zeros(T2,N)

  # Off-Diagonal Components of Khat
  ftmp = copy(f)
  ifresonant = fill(false,n)
  for j in 1:n
    if abs(λc[j] - λe) < 1.0e-12
      ifresonant[j] = true
    end
  end

  Γ = ObliqueSubspaceRemoval3!(ftmp,V,W,B,ifresonant,ngs)

  # Build Extended Matrices
  Ne  = N+1

  resonance = false

  f2 = ftmp'*(B.*ftmp)
  fnorm = sqrt(abs(f2[1]))

  if fnorm > AInp.tol

    # Set the proper Arnoldi vector length
    AInp.vlen     = Ne
    DArnOut       = REPStepArn(L,B,σ,V,W,ifresonant,ftmp,λe,SInp,AInp,lbc,rbc)

    # DArnOut     = StepArn(Le,Be,SInp,AInp,lbc,rbc) 
    ii            = argmin(abs.(DArnOut.evals .- λe))
    λfound        = DArnOut.evals[ii]
    θt            = DArnOut.evecs[Ne,ii]

    for i in LinearIndices(vh)
      vh[i]  = DArnOut.evecs[i,ii]./θt     # ensure extended variable == 1.0
    end  

    if (abs(λe - λfound) > 10*AInp.tol)
      println("Mismatched Eigenvalues: $λe, $λfound")
    else
      println("Eigenvalue Found. $λfound")
    end

  end

  # Adjoint Forcing modes
  wh    = zeros(T2,N)
  z     = zeros(T2,1,n)             # Row Vector
  for i in 1:n
    z[i] = -vh'*(B.*W[:,i])
  end

  extmode = ExtendedMode(λe,Γ,vh,wh,z)

  return extmode
end  
#---------------------------------------------------------------------- 













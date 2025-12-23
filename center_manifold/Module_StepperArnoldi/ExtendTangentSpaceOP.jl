# Extending the tangent space
#---------------------------------------------------------------------- 
function ExtLOP(xe::AbstractVector{T1},OP,B::AbstractVector{T3},f::AbstractVector{T2},ω::T2) where {T1,T2,T3<:Number}

  Ne        = length(xe)
  N         = Ne-1
  v         = zeros(T1,Ne)
  x         = view(xe,1:N)

  v[1:N]    = OP(x) .+ B.*f*xe[Ne]
  v[Ne]     = ω*xe[Ne]

  return v
end  
#---------------------------------------------------------------------- 
function ExtendedTangentSpacesOP(OP,B::AbstractVector{T3},λc::AbstractVector{T1},V::AbstractMatrix{T2},W::AbstractMatrix{T2},F::AbstractMatrix{T2},λe::AbstractVector{T1},restricted::Bool,AInp::ArnoldiInput,SInp::StepperInput,lbc::Bool,rbc::Bool) where {T1,T2,T3<:Number}

  ne = length(λe)

  EM = Vector{StepperArnoldi.ExtendedMode}(undef,ne)

  for i in 1:ne
    f     = view(F,:,i)
    AInp.ifeigshift = true
    AInp.eigshift   = λe[i]
    if (restricted)
      EM[i] = ExtendTangentSpaceRestrictedOP(OP,B,λc,V,W,f,λe[i],AInp,SInp,lbc,rbc)
    else
      EM[i] = ExtendTangentSpaceOP(OP,B,λc,V,W,f,λe[i],AInp,SInp,lbc,rbc)
    end
  end  

  EModes = ExtendedModes(EM)

  return EModes
end
#---------------------------------------------------------------------- 
function ExtendTangentSpaceOP(OP,B::AbstractVector{T3},λc::AbstractVector{T1},V::AbstractMatrix{T2},W::AbstractMatrix{T2},f::AbstractVector{T2},λe::T1,AInp::ArnoldiInput,SInp::StepperInput,lbc::Bool,rbc::Bool) where {T1,T2,T3<:Number}

  ngs = 2
  N   = length(B)
  n   = length(λc)

  Γ   = zeros(T2,n)
  vh  = zeros(T2,N)
  fe  = copy(f)

  # Off-Diagonal Components of Khat
  ftmp = copy(f)
  for ig in 1:ngs
    α = zeros(T2,n)
    for j in 1:n
      if abs(λc[j] - λe) < 1.0e-12
        α[j]  = W[:,j]'*(B.*ftmp)
        Γ[j]  = Γ[j] + α[j]
      end
    end
    ftmp .= ftmp .- V*α
  end

  # Build Extended Matrices
  Ne  = N+1
  Ve  = zeros(T2,Ne,n)
  We  = zeros(T2,Ne,n)
  # fe  = zeros(T2,Ne)

  resonance = false
  nr        = 0
  for j in 1:n
    if abs(λc[j] - λe) < AInp.tol
      # println("Resonant λ: $j")
      resonance = true
      nr  = nr + 1

      di = (nr-1)*Ne + 1    # destination index
      si = (j-1)*N + 1      # source index
      copyto!(Ve,di,V,si,N)
      copyto!(We,di,W,si,N)
    end
  end  
  # @views SEM1D.SEM_SetBC!(ftmp,Inp.lbc,Inp.rbc)

  if norm(ftmp) > AInp.tol
    # Extended Mass
    Be   = ones(T3,Ne)
    copyto!(Be,1,B,1,N)

    # Extended Operator
    # if typeof(L) <: SparseMatrixCSC
    #   Le = spzeros(T2,Ne,Ne)
    # else
    #   Le = zeros(T2,Ne,Ne)
    # end
    # Le[1:N,1:N]   = L
    # Le[1:N,N+1]   = B.*ftmp
    # Le[N+1,N+1]   = λe
    OPLe(x)       = ExtLOP(x,OP,B,ftmp,λe) 

    ifconv = false
    # Set the proper Arnoldi vector length
    AInp.vlen     = Ne
    if (resonance)
      Ve_view     = view(Ve,:,1:nr)
      We_view     = view(We,:,1:nr)
      DArnOut     = RestrictedStepArnOP(OPLe,Be,Ve_view,We_view,SInp,AInp,lbc,rbc) 
    else
      DArnOut     = StepArnOP(OPLe,Be,SInp,AInp,lbc,rbc) 
    end
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
  end       # norm(ftmp)>AInp.tol

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
function ExtendTangentSpaceRestrictedOP(OP,B::AbstractVector{T3},λc::AbstractVector{T1},V::AbstractMatrix{T2},W::AbstractMatrix{T2},f::AbstractVector{T2},λe::T1,AInp::ArnoldiInput,SInp::StepperInput,lbc::Bool,rbc::Bool) where {T1,T2,T3<:Number}

  ngs = 2
  N   = length(B)
  n   = length(λc)

  Γ   = zeros(T2,n)
  vh  = zeros(T2,N)

  # Off-Diagonal Components of Khat
  ftmp = copy(f)
  for ig in 1:ngs
    α = zeros(T2,n)
    for j in 1:n
      α[j]  = W[:,j]'*(B.*ftmp)
      Γ[j]  = Γ[j] + α[j]
    end
    ftmp .= ftmp .- V*α
  end

  # Build Extended Matrices
  Ne  = N+1

  if norm(ftmp) > AInp.tol
    # Extended Mass
    Be   = ones(T3,Ne)
    copyto!(Be,1,B,1,N)

    # Extended Operator
    if typeof(L) <: SparseMatrixCSC
      Le = spzeros(T2,Ne,Ne)
    else
      Le = zeros(T2,Ne,Ne)
    end
    # Le[1:N,1:N]   = L
    # Le[1:N,N+1]   = B.*ftmp
    # Le[N+1,N+1]   = λe
    OPLe(x)       = ExtLOP(x,OP,B,ftmp,λe) 

    Ve  = zeros(T2,Ne,n)
    We  = zeros(T2,Ne,n)
    for j in 1:n
      di = (j-1)*Ne + 1    # destination index
      si = (j-1)*N + 1     # source index
      copyto!(Ve,di,V,si,N)
      copyto!(We,di,W,si,N)
    end  

    ifconv = false
    # Set the proper Arnoldi vector length
    AInp.vlen     = Ne
    DArnOut       = RestrictedStepArn(Le,Be,Ve,We,SInp,AInp,lbc,rbc) 
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













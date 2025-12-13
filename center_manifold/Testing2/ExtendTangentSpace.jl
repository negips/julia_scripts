# Extending the tangent space
#---------------------------------------------------------------------- 
function ExtendTangentSpace(λc::AbstractVector{T1},L::AbstractMatrix{T2},B::AbstractVector{T3},V::AbstractMatrix{T2},W::AbstractMatrix{T2},F::AbstractMatrix{T4},λe::AbstractVector{T1},ArnInp,StpInp,Inp) where {T1,T2,T3,T4<:Number}

  ne = length(λe)

  EM = Vector{StepperArnoldi.ExtendedMode}(undef,ne)

  for i in 1:ne
    f     = view(F,:,i)
    EM[i] = ExtendTangentSpace(λc,L,B,V,W,f,λe[i],ArnInp,StpInp,Inp)
  end  

  EModes = ExtendedModes(EM)

  return EModes
end
#---------------------------------------------------------------------- 
function ExtendTangentSpace(λc::AbstractVector{T1},L::AbstractMatrix{T2},B::AbstractVector{T3},V::AbstractMatrix{T2},W::AbstractMatrix{T2},f::AbstractVector{T4},λe::T1,ArnInp,StpInp,Inp) where {T1,T2,T3,T4<:Number}

  ngs = 2
  N   = size(L,2)
  n   = length(λc)

  Γ   = zeros(T2,n)
  vh  = zeros(T2,N)

  # Off-Diagonal Components of Khat
  ftmp = copy(f)
  for ig in 1:ngs
    α = zeros(T2,n)
    for j in 1:n
      if abs(λc[j] - λe) < 1.0e-12
        α[j]  = W[:,j]'*(B.*f)
        Γ[j]  = Γ[j] + W[:,j]'*(B.*f)
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
    if abs(λc[j] - λe) < ArnInp.tol
      println("Resonant λ: $j")
      resonance = true
      nr  = nr + 1
 
      di = (nr-1)*Ne + 1    # destination index
      si = (j-1)*N + 1      # source index
      copyto!(Ve,di,V,si,N)
      copyto!(We,di,W,si,N)
    end
  end  
  # @views SEM1D.SEM_SetBC!(ftmp,Inp.lbc,Inp.rbc)

  if norm(ftmp) > EArnInp.tol
    # Extended Mass
    Be   = ones(T3,Ne)
    copyto!(Be,1,B,1,N)

    # Extended Operator
    if typeof(L) <: SparseMatrixCSC
      Le = spzeros(T2,Ne,Ne)
    else
      Le = zeros(T2,Ne,Ne)
    end
    Le[1:N,1:N]   = L
    Le[1:N,N+1]   = B.*ftmp
    Le[N+1,N+1]   = λe

    if (resonance)
      Ve_view     = view(Ve,:,1:nr)
      We_view     = view(We,:,1:nr)
      EArnDir     = StepperArnoldi.RestrictedStepArn(Le,Be,Ve_view,We_view,EStpInp,EArnInp,Inp.lbc,Inp.rbc) 
    else
      EArnDir     = StepperArnoldi.StepArn(Le,Be,EStpInp,EArnInp,Inp.lbc,Inp.rbc) 
    end
    ii            = argmin(abs.(EArnDir.evals .- λe))
    θt            = EArnDir.evecs[ndof+1,ii]
    vh           .= EArnDir.evecs[1:ndof,ii]./θt     # ensure extended variable == 1.0
  end

  # Adjoint Forcing modes
  wh    = zeros(T2,N)
  z     = zeros(T2,1,n)             # Row Vector
  for i in 1:n
    z[i] = -vh'*(B.*W[:,i])
  end

  extmode = StepperArnoldi.ExtendedMode(λe,Γ,vh,wh,z)

  return extmode
end  














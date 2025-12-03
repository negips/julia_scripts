# Functions needed for Ginzburg-Landau Center-manifold Evaluation
#----------------------------------------------------------------------  
function Get_SEM1D_Input()

  Np    = 12
  Npd   = 19
  nel   = 61
  xs    = 0.0
  xe    = 40.0
  iflbc = true
  ifrbc = false
  
  # Input parameters
  Inp   = SEM1D.SEMInput(Np,Npd,nel,xs,xe,iflbc,ifrbc)

  return Inp
end  
#---------------------------------------------------------------------- 
function renormalize_evec!(v::AbstractVector{T},j0::Int) where {T}

  vj        = v[j0]
  absv      = abs(vj)
  th        = atan(imag(vj)/absv,real(vj)/absv)
  ph        = π/4.0 - th
  for i in eachindex(v)
    v[i]    = v[i]*exp(ph*im)
  end  

  return nothing
end  
#----------------------------------------------------------------------
function renormalize_evecs!(v::AbstractVector{T},w::AbstractVector{T},B::AbstractVector{S}) where {T<:Number,S<:Number}

  β   = sqrt(v'*(B.*v))
  for i in eachindex(v)
    v[i] = v[i]/β
  end

  β   = w'*(B.*v)
  for i in eachindex(w)
    w[i] = w[i]/(β')
  end

  return nothing
end  
#---------------------------------------------------------------------- 
function ForcingShape(B::AbstractVector{T},xg::AbstractVector{Float64},x0::Float64,σ::Float64) where {T <: Number}

  # Forcing Shape
  ψ     = exp.(-((xg .- x0)/σ).^2)
  ψn    = sqrt(ψ'*(Bg.*ψ))
  ψ    .= ψ./ψn

  return ψ
end
#---------------------------------------------------------------------- 
function SetForcingShape!(ψ::AbstractVector{T},B::AbstractVector{S},xg::AbstractVector{Float64},x0::Float64,σ::Float64) where {T,S <: Number}

  # Forcing Shape
  for i in LinearIndices(ψ)
    ψ[i] = exp(-((xg[i] - x0)/σ)^2)
  end  
  ψn    = sqrt(ψ'*(Bg.*ψ))
  ψ    .= ψ./ψn

  return nothing
end
#---------------------------------------------------------------------- 


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
function Set_GL_Params()
  # δ1  = -1.0 + 0.0im            # -U
  # δ2  =  0.741 + 1.025im        #  μ0
  # δ3  = -0.125 + 0.0im          #  μx
  # δ4  = (1.0 - 1.0im)/sqrt(2.0) #  γ
  # δ5  = (-0.1 + 0.1im)          #  Nonlinear Coefficient 

  # δ     = ones(ComplexF64,5)    #  Parameters
  δ1  = -0.1 + 0.0im            # -U
  # δ2  =  0.741 + 1.025im      #  μ0
  δ3  = -0.125 + 0.0im          #  μx
  δ4  = (1.0 - 0.00im)/sqrt(1.0)    #  γ
  δ5  = (-0.1 + 0.1im)          #  Nonlinear Coefficient 

  ω1  =  0.0 + 1.0im
  #δ   = SEM1D.SetGLParams(ω1,δ5)
  δ   = SEM1D.SetGLParams(ω1,δ1,δ3,δ4,δ5)

  return δ
end  
#---------------------------------------------------------------------- 
function Set_StepperParams()

  ifadjoint         = false
  ifoptimal         = false
  ifverbose         = false
  verbosestep       = 2000
  nsteps            = 2000
  dt                = 2.5e-5
  StpInp            = StepperArnoldi.StepperInput(ifadjoint,ifoptimal,ifverbose,verbosestep,nsteps,dt)

  return StpInp
end  
#---------------------------------------------------------------------- 
function Set_ArnoldiParams()

  ifarnoldi         = true 
  ifverbose         = false
  vlen              = ndof
  nev               = 1
  ekryl             = 15  
  lkryl             = nev + ekryl 
  ngs               = 2
  bsize             = 1
  outer_iterations  = 100
  tol               = 1.0e-12
  ArnInp            = StepperArnoldi.ArnoldiInput(ifarnoldi,ifverbose,vlen,nev,ekryl,lkryl,ngs,bsize,outer_iterations,tol)

  return ArnInp
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
function ForcingParams()

  x0 = 5.0
  κ  = 0.0

  return x0,κ
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
function SetForcingShape!(ψ::AbstractVector{T},B::AbstractVector{S},xg::AbstractVector{Float64},x0::Float64,σ::Float64,κ::Float64) where {T,S <: Number}

  # Forcing Shape
  for i in LinearIndices(ψ)
    xi   = xg[i]
    ψ[i] = exp(-((xi - x0)/σ)^2)*(exp(im*κ*(xi - x0)))
  end  
  ψn    = sqrt(ψ'*(Bg.*ψ))
  ψ    .= ψ./ψn

  return nothing
end
#---------------------------------------------------------------------- 
function Get_AsymptoticFieldx(ind::Int,z::AbstractVector{T},Y1::AbstractMatrix{T},Y2::AbstractMatrix{T},Y3::AbstractMatrix{T}) where {T <: Number}

  val = T(0)

  z1  = z
  val = val + transpose(Y1[ind,:])*z1 

  Ord = 2
  z2  = CenterManifold.EvaluateNonLinear(z,Ord)
  val = val + transpose(Y2[ind,:])*z2

  Ord = 3
  z3  = CenterManifold.EvaluateNonLinear(z,Ord)
  val = val + transpose(Y3[ind,:])*z3

  return val
end
#---------------------------------------------------------------------- 
function Get_AsymptoticField(z::AbstractVector{T},Y1::AbstractMatrix{T},Y2::AbstractMatrix{T},Y3::AbstractMatrix{T}) where {T <: Number}

  ndof,nv = size(Y1)

  fld = zeros(T,ndof)

  z1   = z
  fld .= fld .+ Y1*z1 

  Ord = 2
  z2  = CenterManifold.EvaluateNonLinear(z,Ord)
  fld .= fld .+ Y2*z2

  Ord = 3
  z3  = CenterManifold.EvaluateNonLinear(z,Ord)
  fld .= fld .+ Y3*z3

  return fld
end
#---------------------------------------------------------------------- 
struct Graphics

  lafs::Int       # Label Font Size
  lgfs::Int       # Legend Font Size
  mksz::Int       # Marker Size
  figsz1::Vector{Float64}
  figsz2::Vector{Float64}
  figsz3::Vector{Float64}

end  
#---------------------------------------------------------------------- 
function setgraphics(screen)

  if screen == 1
    # hp spectre
    lafs      = 16        # Label font size
    lgfs      = 12        # Legend font size
    mksz      = 6
    figsz1    = [ 8.0, 6.0]
    figsz2    = [10.0, 7.0]
    figsz3    = [12.0, 5.0]
  elseif screen == 2
    # workstation
    lafs      = 24        # Label font size
    lgfs      = 20        # Legend font size
    mksz      = 12
    figsz1    = [16.0, 14.0]
    figsz2    = [20.0, 16.0]
    figsz3    = [24.0, 12.0]
  end  

  grh = Graphics(lafs,lgfs,mksz,figsz1,figsz2,figsz3)

  return grh 
end
#---------------------------------------------------------------------- 




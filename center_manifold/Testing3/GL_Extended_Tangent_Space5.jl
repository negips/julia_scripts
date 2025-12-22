# Functions for Extending the tangent space
#---------------------------------------------------------------------- 
function GLExtendTangentSpace(L::AbstractMatrix{T1},LC::AbstractMatrix{T1},B::AbstractVector{T2},λsys::AbstractVector{T3},V::AbstractMatrix{T1},W::AbstractMatrix{T1},λe::AbstractVector{T3},Le::AbstractMatrix{T1},restricted::Bool,lbc::Bool,rbc::Bool) where {T1,T2,T3<:Number}

  Nby2      = size(L,2)
  N         = 2*Nby2
  ind1      = 1:Nby2
  ind2      = Nby2+1:N

  nsys      = length(λsys)
  ne        = length(λe)
  m         = nsys + ne

  # Stepper-Arnoldi (for Extended Operator Calculation)
  #-------------------------------------------------- 
  ifadjoint         = false
  ifoptimal         = false
  ifverbose         = false
  verbosestep       = 500
  nsteps            = 500
  dt                = 1.0e-4
  EStpInp           = StepperArnoldi.StepperInput(ifadjoint,ifoptimal,ifverbose,verbosestep,nsteps,dt)
  
  ifarnoldi         = true 
  ifverbose         = false
  ifeigshift        = true
  vlen              = Nby2+1
  nev               = 1
  ekryl             = 15  
  lkryl             = nev + ekryl 
  eigshift          = 0.0 + 0.0im
  ngs               = 2
  bsize             = 1
  outer_iterations  = 100
  tol               = 1.0e-12
  EArnInp           = StepperArnoldi.ArnoldiInput(ifarnoldi,ifverbose,ifeigshift,vlen,nev,ekryl,lkryl,eigshift,ngs,bsize,outer_iterations,tol)
  

  Emodes1    = StepperArnoldi.ExtendedTangentSpaces(L ,B,λsys,V[ind1,:],W[ind1,:],Le[ind1,:],λe,restricted,EArnInp,EStpInp,lbc,rbc);
  Emodes2    = StepperArnoldi.ExtendedTangentSpaces(LC,B,λsys,V[ind2,:],W[ind2,:],Le[ind2,:],λe,restricted,EArnInp,EStpInp,lbc,rbc);

  Z          = Emodes1.Z .+ Emodes2.Z
  Γ          = Emodes1.Γ .+ Emodes2.Γ

  Ve         = zeros(T1,N,ne)
  Ve[ind1,:] = copy(Emodes1.Ve)
  Ve[ind2,:] = copy(Emodes2.Ve)
  
  We         = zeros(T1,N,ne)
  We[ind1,:] = copy(Emodes1.We)
  We[ind2,:] = copy(Emodes2.We)

  Emodes = StepperArnoldi.ExtendedModes(λe,Γ,Ve,We,Z)

  return Emodes
end
#---------------------------------------------------------------------- 
function GLExtendPertTangentSpace(L::AbstractMatrix{T1},LC::AbstractMatrix{T1},B::AbstractVector{T2},λsys::AbstractVector{T3},σ::AbstractVector{T3},V::AbstractMatrix{T1},W::AbstractMatrix{T1},λe::AbstractVector{T3},Le::AbstractMatrix{T1},restricted::Bool,lbc::Bool,rbc::Bool) where {T1,T2,T3<:Number}

  Nby2      = size(L,2)
  N         = 2*Nby2
  ind1      = 1:Nby2
  ind2      = Nby2+1:N

  nsys      = length(λsys)
  ne        = length(λe)
  m         = nsys + ne

  # Stepper-Arnoldi (for Extended Operator Calculation)
  #-------------------------------------------------- 
  ifadjoint         = false
  ifoptimal         = false
  ifverbose         = false
  verbosestep       = 500
  nsteps            = 500
  dt                = 1.0e-4
  EStpInp           = StepperArnoldi.StepperInput(ifadjoint,ifoptimal,ifverbose,verbosestep,nsteps,dt)
  
  ifarnoldi         = true 
  ifverbose         = false
  ifeigshift        = true
  vlen              = Nby2+1
  nev               = 1
  ekryl             = 15  
  lkryl             = nev + ekryl 
  eigshift          = 0.0 + 0.0im
  ngs               = 2
  bsize             = 1
  outer_iterations  = 100
  tol               = 1.0e-12
  EArnInp           = StepperArnoldi.ArnoldiInput(ifarnoldi,ifverbose,ifeigshift,vlen,nev,ekryl,lkryl,eigshift,ngs,bsize,outer_iterations,tol)
  

  Emodes1    = StepperArnoldi.REPTangentSpaces(L ,B,λsys,σ,V[ind1,:],W[ind1,:],Le[ind1,:],λe,restricted,EArnInp,EStpInp,lbc,rbc);
  Emodes2    = StepperArnoldi.REPTangentSpaces(LC,B,λsys,σ,V[ind2,:],W[ind2,:],Le[ind2,:],λe,restricted,EArnInp,EStpInp,lbc,rbc);

  Z          = Emodes1.Z .+ Emodes2.Z
  Γ          = Emodes1.Γ .+ Emodes2.Γ

  Ve         = zeros(T1,N,ne)
  Ve[ind1,:] = copy(Emodes1.Ve)
  Ve[ind2,:] = copy(Emodes2.Ve)
  
  We         = zeros(T1,N,ne)
  We[ind1,:] = copy(Emodes1.We)
  We[ind2,:] = copy(Emodes2.We)

  Emodes = StepperArnoldi.ExtendedModes(λe,Γ,Ve,We,Z)

  return Emodes
end
#---------------------------------------------------------------------- 














# Functions for Extending the tangent space
#---------------------------------------------------------------------- 
function GLExtendTangentSpace2(L::AbstractMatrix{T1},LC::AbstractMatrix{T1},B::AbstractVector{T2},λc::AbstractVector{T3},V::AbstractMatrix{T1},W::AbstractMatrix{T1},λp::AbstractVector{T3},Lp::AbstractMatrix{T1},λh::AbstractVector{T3},Lh::AbstractMatrix{T1},restricted::Bool,lbc::Bool,rbc::Bool) where {T1,T2,T3<:Number}

  Nby2            = size(L,2)
  N               = 2*Nby2
  ind1            = 1:Nby2
  ind2            = Nby2+1:N

  n               = length(λc)
  p               = length(λp)
  h               = length(λh)

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
  

  PEmodes1   = StepperArnoldi.ExtendedTangentSpaces(L ,B,λc,V[ind1,:],W[ind1,:],Lp[ind1,:],λp,restricted,EArnInp,EStpInp,lbc,rbc);
  PEmodes2   = StepperArnoldi.ExtendedTangentSpaces(LC,B,λc,V[ind2,:],W[ind2,:],Lp[ind2,:],λp,restricted,EArnInp,EStpInp,lbc,rbc);
  HEmodes1   = StepperArnoldi.ExtendedTangentSpaces(L ,B,λc,V[ind1,:],W[ind1,:],Lh[ind1,:],λh,restricted,EArnInp,EStpInp,lbc,rbc);
  HEmodes2   = StepperArnoldi.ExtendedTangentSpaces(LC,B,λc,V[ind2,:],W[ind2,:],Lh[ind2,:],λh,restricted,EArnInp,EStpInp,lbc,rbc);

  Zp         = PEmodes1.Z .+ PEmodes2.Z
  Zh         = HEmodes1.Z .+ HEmodes2.Z
  Γp         = PEmodes1.Γ .+ PEmodes2.Γ
  Γh         = HEmodes1.Γ .+ HEmodes2.Γ

  Vp         = zeros(T1,N,p)
  Vp[ind1,:] = copy(PEmodes1.Ve)
  Vp[ind2,:] = copy(PEmodes2.Ve)

  Vh         = zeros(T1,N,h)
  Vh[ind1,:] = copy(HEmodes1.Ve)
  Vh[ind2,:] = copy(HEmodes2.Ve)
  
  Wp         = zeros(T1,N,p)
  Wp[ind1,:] = copy(PEmodes1.We)
  Wp[ind2,:] = copy(PEmodes2.We)

  Wh         = zeros(T1,N,h)
  Wh[ind1,:] = copy(HEmodes1.We)
  Wh[ind2,:] = copy(HEmodes2.We)

  PEmodes = StepperArnoldi.ExtendedModes(λp,Γp,Vp,Wp,Zp)
  HEmodes = StepperArnoldi.ExtendedModes(λh,Γh,Vh,Wh,Zh)

  #return PEmodes1,PEmodes2,HEmodes1,HEmodes2
  return PEmodes,HEmodes
end
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
# function GLExtendTangentSpace2OP(OP,OPC,B::AbstractVector{T2},λc::AbstractVector{T3},V::AbstractMatrix{T1},W::AbstractMatrix{T1},λe::AbstractVector{T3},Le::AbstractMatrix{T1},restricted::Bool,lbc::Bool,rbc::Bool) where {T1,T2,T3<:Number}
# 
#   Nby2      = length(B)
#   N         = 2*Nby2
#   ind1      = 1:Nby2
#   ind2      = Nby2+1:N
# 
#   n         = length(λc)
#   ne        = length(λe)
# 
#   # Stepper-Arnoldi (for Extended Operator Calculation)
#   #-------------------------------------------------- 
#   ifadjoint         = false
#   ifoptimal         = false
#   ifverbose         = false
#   verbosestep       = 500
#   nsteps            = 500
#   dt                = 1.0e-4
#   EStpInp           = StepperArnoldi.StepperInput(ifadjoint,ifoptimal,ifverbose,verbosestep,nsteps,dt)
#   
#   ifarnoldi         = true 
#   ifverbose         = false
#   ifeigshift        = true
#   vlen              = Nby2+1
#   nev               = 1
#   ekryl             = 15  
#   lkryl             = nev + ekryl 
#   eigshift          = 0.0 + 0.0im
#   ngs               = 2
#   bsize             = 1
#   outer_iterations  = 100
#   tol               = 1.0e-12
#   EArnInp           = StepperArnoldi.ArnoldiInput(ifarnoldi,ifverbose,ifeigshift,vlen,nev,ekryl,lkryl,eigshift,ngs,bsize,outer_iterations,tol)
#   
# 
#   Emodes1    = StepperArnoldi.ExtendedTangentSpacesOP(OP,B,λc,V[ind1,:],W[ind1,:],Le[ind1,:],λe,restricted,EArnInp,EStpInp,lbc,rbc);
#   Emodes2    = StepperArnoldi.ExtendedTangentSpacesOP(OPC,B,λc,V[ind2,:],W[ind2,:],Le[ind2,:],λe,restricted,EArnInp,EStpInp,lbc,rbc);
# 
#   Z          = Emodes1.Z .+ Emodes2.Z
#   Γ          = Emodes1.Γ .+ Emodes2.Γ
# 
#   Ve         = zeros(T1,N,ne)
#   Ve[ind1,:] = copy(Emodes1.Ve)
#   Ve[ind2,:] = copy(Emodes2.Ve)
#   
#   We         = zeros(T1,N,ne)
#   We[ind1,:] = copy(Emodes1.We)
#   We[ind2,:] = copy(Emodes2.We)
# 
#   Emodes = StepperArnoldi.ExtendedModes(λe,Γ,Ve,We,Z)
# 
#   return Emodes
# end
#---------------------------------------------------------------------- 
function GLExtendTangentSpace2OP(OP,OPC,B::AbstractVector{T2},λc::AbstractVector{T3},V::AbstractMatrix{T1},W::AbstractMatrix{T1},λe::AbstractVector{T3},Le::AbstractMatrix{T1},restricted::Bool,lbc::Bool,rbc::Bool) where {T1,T2,T3<:Number}

  Nby2      = length(B)
  N         = 2*Nby2
  ind1      = 1:Nby2
  ind2      = Nby2+1:N

  n         = length(λc)
  ne        = length(λe)

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
  

  Emodes1    = StepperArnoldi.ExtendedTangentSpacesOP(OP,B,λc,V[ind1,:],W[ind1,:],Le[ind1,:],λe,restricted,EArnInp,EStpInp,lbc,rbc);
  Emodes2    = StepperArnoldi.ExtendedTangentSpacesOP(OPC,B,λc,V[ind2,:],W[ind2,:],Le[ind2,:],λe,restricted,EArnInp,EStpInp,lbc,rbc);

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
function GetExternalForcing(x::AbstractVector{T1},B::AbstractVector{T2},lbc::Bool,rbc::Bool) where {T1,T2<:Number}

  # Forcing Shape
  x0,κ  = ForcingParams()
  Nby2  = length(B)
  N     = 2*Nby2
  ψ     = zeros(ComplexF64,Nby2)
  σ     = 1.0
  SetForcingShape!(ψ,B,x,x0,σ,κ)
  SEM1D.SEM_SetBC!(ψ,lbc,rbc)

  h     = 2
  Lθ    = zeros(ComplexF64,N,h)
  f2    = [ψ;  vzro]
  f3    = [vzro;ψ]
  
  Lθ    = [ψ                        zeros(ComplexF64,Nby2);
           zeros(ComplexF64,Nby2)   conj.(ψ)]

  if ifresonant
    λh    = [1.0im; -1.0im;]
  else
    λh    = [1.3im; -1.3im;]
  end


  return ψ,Lθ,λh
end
#----------------------------------------------------------------------














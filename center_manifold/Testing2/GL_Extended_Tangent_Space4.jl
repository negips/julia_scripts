# Testing the module

# Extending the tangent space
#---------------------------------------------------------------------- 
println("Extended Tangent Space using Arnoldi.")

function GLExtendTangentSpace2(L::AbstractMatrix{T1},LC::AbstractMatrix{T1},B::AbstractVector{T2},λc::AbstractVector{T3},V::AbstractMatrix{T1},W::AbstractMatrix{T1},λp::AbstractVector{T3},Lp::AbstractMatrix{T1},λh::AbstractVector{T3},Lh::AbstractMatrix{T1},lbc::Bool,rbc::Bool) where {T1,T2,T3<:Number}

  Nby2      = size(L,2)
  N         = 2*Nby2
  ind1      = 1:Nby2
  ind2      = Nby2+1:N

  n         = length(λc)
  p         = length(λp)
  h         = length(λh)

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
  

  PEmodes1   = StepperArnoldi.ExtendedTangentSpaces(L ,B,λc,V[ind1,:],W[ind1,:],Lp[ind1,:],λp,EArnInp,EStpInp,lbc,rbc);
  PEmodes2   = StepperArnoldi.ExtendedTangentSpaces(LC,B,λc,V[ind2,:],W[ind2,:],Lp[ind2,:],λp,EArnInp,EStpInp,lbc,rbc);
  HEmodes1   = StepperArnoldi.ExtendedTangentSpaces(L ,B,λc,V[ind1,:],W[ind1,:],Lh[ind1,:],λh,EArnInp,EStpInp,lbc,rbc);
  HEmodes2   = StepperArnoldi.ExtendedTangentSpaces(LC,B,λc,V[ind2,:],W[ind2,:],Lh[ind2,:],λh,EArnInp,EStpInp,lbc,rbc);

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
function GLExtendTangentSpace2(L::AbstractMatrix{T1},LC::AbstractMatrix{T1},B::AbstractVector{T2},λc::AbstractVector{T3},V::AbstractMatrix{T1},W::AbstractMatrix{T1},λe::AbstractVector{T3},Le::AbstractMatrix{T1},lbc::Bool,rbc::Bool) where {T1,T2,T3<:Number}

  Nby2      = size(L,2)
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
  

  Emodes1    = StepperArnoldi.ExtendedTangentSpaces(L ,B,λc,V[ind1,:],W[ind1,:],Le[ind1,:],λe,EArnInp,EStpInp,lbc,rbc);
  Emodes2    = StepperArnoldi.ExtendedTangentSpaces(LC,B,λc,V[ind2,:],W[ind2,:],Le[ind2,:],λe,EArnInp,EStpInp,lbc,rbc);

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

ifresonant = false
emodeplot  = true

Nby2  = ArnInp.vlen
N     = Nby2*2
n     = length(λc)
p     = 2
h     = 2
m     = n+p+h
ind1  = 1:Nby2
ind2  = Nby2+1:N

Bg2   = [Bg; Bg]
Bg2M  = diagm(Bg2)

Lν    = zeros(ComplexF64,N,p)
λν    = zeros(ComplexF64,p)
# Forcing Shape
x0,κ  = ForcingParams()
ψ     = zeros(ComplexF64,Nby2)
σ     = 1.0
SetForcingShape!(ψ,Bg,xg,x0,σ,κ)
SEM1D.SEM_SetBC!(ψ,Inp.lbc,Inp.rbc)
Lθ    = zeros(ComplexF64,N,h)
f1    = 1.0/(sqrt(2.0))*[ψ;ψ]
f2    = [ψ;  vzro]
f3    = [vzro;ψ]
Lθ    = [f2 f3]
if ifresonant
  λh    = [1.0im; -1.0im;]
else
  λh    = [1.2im; -1.2im;]
end
# Forcing Modes
ax2.plot(xg,real.(ψ) ,linewidth=2,linestyle="-", color=cm(2),label=L"\mathfrak{R}(ψ)")
ax2.plot(xg,imag.(ψ) ,linewidth=2,linestyle="--",color=cm(2),label=L"\mathfrak{Im}(ψ)")

Λh    = diagm(λh)

# PE,HE = GLExtendTangentSpace2(OPg,OPCg,Bg,λc,V,W,λν,Lν,λh,Lθ,Inp.lbc,Inp.rbc)

λext = [λν[:]; λh[:]]
Lext = zeros(ComplexF64,N,p+h)
for i in 1:p
  for j in 1:N
    Lext[j,i] = Lν[j,i]
  end
end
for i in 1:h
  for j in 1:N
    Lext[j,i+p] = Lθ[j,i]
  end
end
EM = GLExtendTangentSpace2(OPg,OPCg,Bg,λc,V,W,λext,Lext,Inp.lbc,Inp.rbc)

for i in 1:length(λext)
  vnorm = norm(EM.Ve[ind1,i])
  # Plot Mode
  if (emodeplot) && vnorm > 0.0
    j = n+i+20
    ax2.plot(xg,real.(EM.Ve[ind1,i]),linewidth=2,linestyle="-", color=cm(j-1),label=L"\mathfrak{R}(ϕ_{%$j})")
    ax2.plot(xg,imag.(EM.Ve[ind1,i]),linewidth=2,linestyle="--",color=cm(j-1),label=L"\mathfrak{Im}(ϕ_{%$j})")
  end
end  

Vext        = [V EM.Ve]
Wext        = [W EM.We]
Γe          = EM.Γ
Ze          = EM.Z
Λe          = diagm(EM.λe)
Λc          = diagm(λc)
Zero_ne_n   = zeros(ComplexF64,p+h,n)

Khat = [Λc        Γe;
        Zero_ne_n Λe]

Vhat  = [V              EM.Ve;
         Zero_ne_n      I]
What  = [Wext;
         Ze      I]

if (emodeplot)
  ax2.legend(ncols=4,fontsize=Grh.lgfs)
else  
  ax2.legend(ncols=3,fontsize=Grh.lgfs)
end  


# Extended Adjoint Tangent Space
#-------------------------------------------------- 

# println("Extended Tangent Space Done.")


println("Extended Tangent Space (Arnoldi) Done.")














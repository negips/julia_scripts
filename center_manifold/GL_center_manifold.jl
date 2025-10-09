# Testing the module

include("Module_SEM1D/SEM1D.jl")
#using .SEM1D

include("Module_StepperArnoldi/StepperArnoldi.jl")
#using .StepperArnoldi

using LinearAlgebra
using SparseArrays
using Printf
using PolynomialBases


using PyPlot


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
    w[i] = w[i]/β
  end

  return nothing
end  
#---------------------------------------------------------------------- 

Np    = 12
Npd   = 19
nel   = 41
xs    = 0.0
xe    = 40.0
lbc   = true
rbc   = false

# Input parameters
Inp   = SEM1D.SEMInput(Np,Npd,nel,xs,xe,lbc,rbc)
# Nodal Bases
B0    = LobattoLegendre(Inp.N)
Bd    = LobattoLegendre(Inp.Nd)
# Geometric Matrices
GeoM  = SEM1D.SEMGeoMat(B0,Bd,Inp)

δ     = ones(ComplexF64,4)    # Parameters
δ[1]  = -1.0                  # -U
δ[2]  =  0.741 + 1.025im      #  μ0
δ[3]  = -0.125                #  μx
δ[4]  = (1.0 - im)/sqrt(2.0)  #  γ
ifconj = false
if (ifconj)
  δ = con.(δ)
end  

# GinzburgLandau Linear Operators
L,  B, OP,  Conv,   Src, Lap  = SEM1D.GinzburgLandauSparse(δ,Inp,GeoM,B0)
AL, B, AOP, AConv, ASrc, ALap = SEM1D.AdjointGinzburgLandauSparse(δ,Inp,GeoM,B0)

ifperiodic  = false
ndof, glnum = SEM1D.SEM_Global_Num(GeoM.xm1,ifperiodic)
Q,QT        = SEM1D.SEM_QQT(glnum)
vmult       = Q*QT*ones(Float64,length(GeoM.xm1[:]))
vimult      = 1.0./vmult

SEM1D.GLSetBC!(L ,Inp.lbc,Inp.rbc,ifperiodic)
SEM1D.GLSetBC!(AL,Inp.lbc,Inp.rbc,ifperiodic)

Lg    = QT*L*Q
ALg   = QT*AL*Q
Bg    = QT*B
Bgi   = 1.0./Bg
OPg   = Bgi.*Lg
AOPg  = Bgi.*ALg

# Stepper-Arnoldi
#-------------------------------------------------- 
ifadjoint         = false
ifoptimal         = false
ifverbose         = false
verbosestep       = 500
nsteps            = 500
dt                = 1.0e-4
StpInp            = StepperArnoldi.StepperInput(ifadjoint,ifoptimal,ifverbose,verbosestep,nsteps,dt)

ifarnoldi         = true 
ifverbose         = true
vlen              = ndof
nev               = 1
ekryl             = 15  
lkryl             = nev + ekryl 
ngs               = 2
bsize             = 1
outer_iterations  = 100
tol               = 1.0e-10
ArnInp            = StepperArnoldi.ArnoldiInput(ifarnoldi,ifverbose,vlen,nev,ekryl,lkryl,ngs,bsize,outer_iterations,tol)

Ω     = SEM1D.GLAnalyticalSpectra(δ)
if (StpInp.ifadjoint)
  Ω   = conj.(Ω)
end  

close("all")
h1    = figure(num=1,figsize=[8.,6.]);
ax1   = gca()
pl1   = ax1.plot(imag.(Ω),real.(Ω),linestyle="none",marker="o",markersize=8)
ax1.set_ylim([-2.75,0.5])
if (StpInp.ifadjoint)
  ax1.set_xlim([-1.80,0.0])
else  
  ax1.set_xlim([0.0,1.80])
end  

ArnDir      = StepperArnoldi.StepArn( OPg,Bg,StpInp,ArnInp,Inp.lbc,Inp.rbc)
ArnAdj      = StepperArnoldi.StepArn(AOPg,Bg,StpInp,ArnInp,Inp.lbc,Inp.rbc)

pl3   = ax1.plot(imag.(ArnDir.evals),real.(ArnDir.evals),linestyle="none",marker="s",markersize=4)

cm    = get_cmap("tab10")
xg    = QT*(vimult.*GeoM.xm1[:])

# Build Center-Manifold matrices
v1    = copy(ArnDir.evecs[:,1])
w1    = copy(ArnAdj.evecs[:,1])
j0    = 76
renormalize_evec!(v1,j0)
renormalize_evecs!(v1,w1,Bg)
v2    = conj.(v1)
w2    = conj.(w1)
λc    = [im; -im;]

# Forcing
ψ     = exp.(-(xg .- 6.0).^2) 
ψn    = sqrt(ψ'*(Bg.*ψ))
ψ    .= ψ./ψn

h2    = figure(num=2,figsize=[8.,6.])
ax2   = gca()
ax2.plot(xg,real.(v1),linewidth=2,linestyle="-", color=cm(0))
ax2.plot(xg,imag.(v1),linewidth=2,linestyle="--",color=cm(0))
ax2.plot(xg,real.(w1),linewidth=2,linestyle="-", color=cm(1))
ax2.plot(xg,imag.(w1),linewidth=2,linestyle="--",color=cm(1))
ax2.plot(xg,ψ,linewidth=2,linestyle="-",color=cm(2))

zro   = 0.0*v1

V     = [v1       zro;
         zro      v2]
W     = [w1       zro;
         zro      w2]

Nby2  = ArnInp.vlen
N     = Nby2*2
n     = 2
p     = 2
h     = 5

# Parameter modes
Lν    = zeros(ComplexF64,N,p)
ΓP    = zeros(ComplexF64,n,p)
λp    = [0.0im; 0.0im]
for i in 1:p
  for j in 1:n
    if abs(λc[j] - λp[i]) < 1.0e-12
      ΓH[j,i]   = W[:,j]'*(Bg2.*Lν[:,i])
    end
  end
end

# Forcing Modes
ΓH    = zeros(ComplexF64,n,h)
λh    = [0.0im; im; -im; 2.3im; -2.3im]
f1    = 1.0/(sqrt(2.0))*[ψ;ψ]
f2    = [ψ;zro]
f3    = [zro;ψ]
f4    = [ψ;zro]
f5    = [zro;ψ]
F     = [f1 f2 f3 f4 f5]

Bg2   = [Bg;Bg]
for i in 1:h
  for j in 1:n
    if abs(λc[j] - λh[i]) < 1.0e-12
      ΓH[j,i]   = W[:,j]'*(Bg2.*F[:,i])
    end
  end
end

# Reduced Matrix
Khat  = [diagm(λc)                  ΓP                      ΓH;
         zeros(ComplexF64,p,n)      diagm(λp)               zeros(ComplexF64,p,h);
         zeros(ComplexF64,h,n)      zeros(ComplexF64,h,p)   diagm(λh)]

# Extended Eigenspace
vp1   = zeros(ComplexF64,Nby2)
vp2   = conj.(vp1)

Vp    = [vp1 zro;
         zro vp2]

Vh    = zeros(ComplexF64,N,h)
Rh    = zeros(ComplexF64,N,h)
for i in 1:h
  r = copy(F[:,i])
  for j in 1:n
    if abs(λc[j] - λh[i]) < 1.0e-12
      # r .= r .- V[:,j]*ΓH[j,i]
    end
  end  
  vh              = copy(r[1:Nby2])
  Rh[1:Nby2,i]    = copy(r[1:Nby2])
  gmres!(vh,Lg,r[1:Nby2])
  Vh[1:Nby2,i]    = copy(vh)
end  






println("Done.")








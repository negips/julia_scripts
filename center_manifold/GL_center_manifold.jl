# Testing the module

include("Module_SEM1D/SEM1D.jl")
#using .SEM1D

include("Module_StepperArnoldi/StepperArnoldi.jl")
#using .StepperArnoldi

include("Module_CenterManifold/CenterManifold.jl")
#using .CenterManifold

using LinearAlgebra
using SparseArrays
using Printf
using PolynomialBases
using IterativeSolvers

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
    w[i] = w[i]/(β')
  end

  return nothing
end  
#---------------------------------------------------------------------- 
function GLStdAsympNLTerm(n0::Int,Ord0::Int,nc::Int,MV1::AbstractMatrix{T},MV2::AbstractMatrix{T},MV3::AbstractMatrix{T},Ord1::Int,Ord2::Int,Ord3::Int) where {T}

  ind0 = CenterManifold.GetPolynomialIndices(n0,Ord0,nc) .+ 1

  Nt1  = CenterManifold.NInteractionTerms(Ord1,nc)
  Nt2  = CenterManifold.NInteractionTerms(Ord2,nc)
  Nt3  = CenterManifold.NInteractionTerms(Ord3,nc)

  N,c  = size(MV1)
  Nby2 = Int(N/2)

  h    = zeros(T,N)

  topi = 1:Nby2
  boti = Nby2+1:N

  if (Ord1 + Ord2 + Ord3 == Ord0) && Ord0>=3

    # The standard nonlinearity δ5*|A*A|A and δ5'|A*A|A* only starts at third-order.
    for i in 1:Nt1
      ind1 = CenterManifold.GetPolynomialIndices(i,Ord1,nc) .+ 1
      for j in 1:Nt2
        ind2 = CenterManifold.GetPolynomialIndices(j,Ord2,nc) .+ 1
        for k in 1:Nt3
          ind3 = CenterManifold.GetPolynomialIndices(k,Ord3,nc) .+ 1

          ind_total = [ind1[:];ind2[:];ind3[:]]

          sort!(ind_total)

          if (ind_total == ind0)    # Matching polynomials
            # println("Index Match: $ind1, $ind2, $ind3: $ind0")
            h[topi] .= h[topi] .+ MV1[boti,i].*MV2[topi,j].*MV3[topi,k]
            h[boti] .= h[boti] .+ MV1[topi,i].*MV2[boti,j].*MV3[boti,k]
          end            
        end       # k
      end         # j
    end           # i
  end             # if (Ord1 + Ord2 + Ord3 == Ord0) && Ord0>=3

  return h
end
#---------------------------------------------------------------------- 
function GLLapNLTerm(n0::Int,Ord0::Int,nc::Int,MV1::AbstractMatrix{T},Ord1::Int,mode::Int,modeC::Int,Lap::AbstractMatrix{S},LapC::AbstractMatrix{S}) where {T,S}

  # mode    - Mode No of δ4
  # modeC   - Mode No of δ4' (conjugate)

  ind0 = CenterManifold.GetPolynomialIndices(n0,Ord0,nc) .+ 1

  Nt1  = CenterManifold.NInteractionTerms(Ord1,nc)

  N,c  = size(MV1)
  Nby2 = Int(N/2)

  h    = zeros(T,N)

  topi = 1:Nby2
  boti = Nby2+1:N

  if (Ord1 + 1 == Ord0) && Ord0>=2

    # The Laplacian nonlinearity {δ4(∇A); δ4'∇(A*)}
    for i in 1:Nt1
      ind1 = CenterManifold.GetPolynomialIndices(i,Ord1,nc) .+ 1

      ind_total = [ind1[:]; mode]
      sort!(ind_total)
      if (ind_total == ind0)    # Matching polynomials
        println("Index Match: $ind1, $mode: $ind0")
        h[topi] .= h[topi] .+ Lap*MV1[topi,i]
      end

      ind_total = [ind1[:]; modeC]
      sort!(ind_total)
      if (ind_total == ind0)    # Matching polynomials
        println("Index Match: $ind1, $modeC: $ind0")
        h[boti] .= h[boti] .+ LapC*MV1[boti,i]
      end            
    end           # i
  end             # if (Ord1 + 1 == Ord0) && Ord0>=2

  return h
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

δ     = ones(ComplexF64,5)    #  Parameters
δ[1]  = -1.0                  # -U
δ[2]  =  0.741 + 1.025im      #  μ0
δ[3]  = -0.125                #  μx
δ[4]  = (1.0 - im)/sqrt(2.0)  #  γ
δ[5]  = (-0.1 + 0.1im)        #  Nonlinear Coefficient 
δc    = conj.(δ)

# GinzburgLandau Linear Operators
L,  B, OP,  Conv,  Src,  Lap  = SEM1D.GinzburgLandauSparse(δ, Inp,GeoM,B0)
LC, B, OPC, ConvC, SrcC, LapC = SEM1D.GinzburgLandauSparse(δc,Inp,GeoM,B0)
AL, B, AOP, AConv, ASrc, ALap = SEM1D.AdjointGinzburgLandauSparse(δ,Inp,GeoM,B0)

ifperiodic  = false
ndof, glnum = SEM1D.SEM_Global_Num(GeoM.xm1,ifperiodic)
Q,QT        = SEM1D.SEM_QQT(glnum)
vmult       = Q*QT*ones(Float64,length(GeoM.xm1[:]))
vimult      = 1.0./vmult

SEM1D.GLSetBC!(L ,Inp.lbc,Inp.rbc,ifperiodic)
SEM1D.GLSetBC!(LC,Inp.lbc,Inp.rbc,ifperiodic)
SEM1D.GLSetBC!(AL,Inp.lbc,Inp.rbc,ifperiodic)

Lg    = QT*L*Q
LCg   = QT*LC*Q

ALg   = QT*AL*Q
Bg    = QT*B
Bgi   = 1.0./Bg
OPg   = Bgi.*Lg
OPCg  = Bgi.*LCg
AOPg  = Bgi.*ALg

Lapg  = Bgi.*(QT*(1.0/δ[4])*Lap*Q)
LapCg = Bgi.*(QT*(1.0/δ[4]')*LapC*Q)


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
tol               = 1.0e-12
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

cm    = get_cmap("tab20")
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

h2    = figure(num=2,figsize=[12.,9.])
ax2   = gca()
ax2.plot(xg,real.(v1),linewidth=2,linestyle="-", color=cm(0),label=L"\mathfrak{R}(ϕ)")
ax2.plot(xg,imag.(v1),linewidth=2,linestyle="--",color=cm(0),label=L"\mathfrak{Im}(ϕ)")
ax2.plot(xg,real.(w1),linewidth=2,linestyle="-", color=cm(1),label=L"\mathfrak{R}(χ)")
ax2.plot(xg,imag.(w1),linewidth=2,linestyle="--",color=cm(1),label=L"\mathfrak{Im}(χ)")
ax2.plot(xg,ψ,linewidth=2,linestyle="-",color=cm(2),label=L"\mathfrak{R}(ψ)")

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
m     = n+p+h

Bg2   = [Bg;Bg]

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
# Parameter modes
Vp    = zeros(ComplexF64,N,p)
Rp    = zeros(ComplexF64,N,p)
ind1  = 1:Nby2
ind2  = Nby2+1:N
for i in 1:p
  r  = copy(Lν[:,i])
  DQ = Matrix{ComplexF64}(I,Nby2,Nby2)
  AQ = Matrix{ComplexF64}(I,Nby2,Nby2)
  for j in 1:n
    if abs(λc[j] - λp[i]) < 1.0e-12
      println("Resonant λp: $i, $j")
      r  .= r .- V[:,j]*ΓP[j,i]
      DQ .= DQ - V[ind1,j]*(W[ind1,j].*Bg)'
      AQ .= AQ - V[ind2,j]*(W[ind2,j].*Bg)'
    end
  end  
  Rp[:,i]         = copy(r)
  
  ω               = λp[i]
  Res1            = DQ*(ω*I - OPg)*DQ
  vp1             = copy(r[ind1])
  gmres!(vp1,Res1,r[ind1])
  Vp[ind1,i]      = copy(vp1)

  Res2            = AQ*(ω*I - OPCg)*AQ
  vp2             = copy(r[ind2])
  gmres!(vp2,Res2,r[ind2])
  Vp[ind2,i]      = copy(vp2)

  vnorm = sqrt(abs(Vp[:,i]'*(Bg2.*Vp[:,i])))
  if vnorm>1.0e-12
    leg = L"\mathfrak{R}(vp"*"$i"*L")"
    ax2.plot(xg,real.(vp1),linewidth=2,linestyle="-", color=cm(n+i),label=leg)
    leg = L"\mathfrak{R}(vp"*"$i"*L")"
    ax2.plot(xg,real.(vp2),linewidth=2,linestyle="--",color=cm(n+i),label=leg)
  end  
end  


# Forcing Modes
Vh    = zeros(ComplexF64,N,h)
Rh    = zeros(ComplexF64,N,h)
ind1  = 1:Nby2
ind2  = Nby2+1:N
for i in 1:h
  r  = copy(F[:,i])
  DQ = Matrix{ComplexF64}(I,Nby2,Nby2)
  AQ = Matrix{ComplexF64}(I,Nby2,Nby2)
  for j in 1:n
    if abs(λc[j] - λh[i]) < 1.0e-12
      println("Resonant λh: $i, $j")
      r  .= r .- V[:,j]*ΓH[j,i]
      DQ .= DQ - V[ind1,j]*(W[ind1,j].*Bg)'
      AQ .= AQ - V[ind2,j]*(W[ind2,j].*Bg)'
    end
  end  
  Rh[:,i]         = copy(r)
  
  ω               = λh[i]
  Res1            = DQ*(ω*I - OPg)*DQ
  vh1             = copy(r[ind1])
  gmres!(vh1,Res1,r[ind1])
  Vh[ind1,i]      = copy(vh1)

  Res2            = AQ*(ω*I - OPCg)*AQ
  vh2             = copy(r[ind2])
  gmres!(vh2,Res2,r[ind2])
  Vh[ind2,i]      = copy(vh2)

  vnorm = sqrt(abs(Vh[:,i]'*(Bg2.*Vh[:,i])))
  if vnorm>1.0e-12
    leg = L"\mathfrak{R}(vh"*"$i"*L")"
    ax2.plot(xg,real.(vh1),linewidth=1,linestyle="-", color=cm(n+p+i),label=leg)
    leg = L"\mathfrak{R}(vh"*"$i"*L")"
    ax2.plot(xg,real.(vh2),linewidth=3,linestyle="--",color=cm(n+p+i),label=leg)
  end  
end  

Vext  = [V  Vp    Vh]

# Second Order Asymptotic terms
Ord   = 2
Nt    = CenterManifold.NInteractionTerms(Ord,m)
H2    = zeros(ComplexF64,N,Nt)

for i in 1:Nt
  ind = CenterManifold.GetPolynomialIndices(i,Ord,m)
  local i1  = ind[1]+1
  local i2  = ind[2]+1

  # println("($i1,$i2)")

  local lap_mode  = n+1
  local lap_cmode = n+2
  Ord1 = 1
  h_asymp   = GLLapNLTerm(i,Ord,m,Vext,Ord1,lap_mode,lap_cmode,Lapg,LapCg)
  # @views SEM1D.SEM_SetBC!(h_asymp[ind1],Inp.lbc,Inp.rbc)
  # @views SEM1D.SEM_SetBC!(h_asymp[ind2],Inp.lbc,Inp.rbc)

  H2[:,i]   = copy(h_asymp)
  hnorm     = h_asymp'*(Bg2.*h_asymp)
  println("Index: $ind, Norm: $hnorm")
  # if abs(hnorm)>1.0e-12
  #   leg = L"\mathfrak{R}(lapl"*"$i"*L")"
  #   ax2.plot(xg,real.(H2[ind1,i]),linewidth=1,linestyle="-", color=cm(m+i),label=leg)
  #   leg = L"\mathfrak{R}(lapl"*"$i"*L")"
  #   ax2.plot(xg,real.(H2[ind2,i]),linewidth=3,linestyle="--",color=cm(m+i),label=leg)
  # end  
end


# Third Order Asymptotic terms
Ord   = 3
Nt    = CenterManifold.NInteractionTerms(Ord,m)
H3    = zeros(ComplexF64,N,Nt)

for i in 1:Nt
  ind = CenterManifold.GetPolynomialIndices(i,Ord,m)
  local i1  = ind[1]+1
  local i2  = ind[2]+1
  local i3  = ind[3]+1

  # println("($i1,$i2)")

  Ord1 = 1
  Ord2 = 1
  Ord3 = 1
  h_asymp   = GLStdAsympNLTerm(i,Ord,m,Vext,Vext,Vext,Ord1,Ord2,Ord3)
  H3[:,i]   = copy(h_asymp)
  hnorm     = h_asymp'*(Bg2.*h_asymp)
  # println("Index: $ind, Norm: $hnorm")
end

ax2.legend()
println("Done.")
















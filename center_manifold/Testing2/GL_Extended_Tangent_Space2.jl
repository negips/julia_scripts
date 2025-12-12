# Testing the module

# Extending the tangent space
#---------------------------------------------------------------------- 

ifresonant = false
emodeplot  = true

Nby2  = ArnInp.vlen
N     = Nby2*2
n     = 2
p     = 2
h     = 2
m     = n+p+h

Bg2   = [Bg; Bg]
Bg2M  = diagm(Bg2)

# Parameter modes
Lν    = zeros(ComplexF64,N,p)
ΓP    = zeros(ComplexF64,n,p)
λp    = zeros(ComplexF64,p)
Λp    = diagm(λp)
for i in 1:p
  for j in 1:n
    if abs(λc[j] - λp[i]) < 1.0e-12
      ΓH[j,i]   = W[:,j]'*(Bg2.*Lν[:,i])
    end
  end
end

# Forcing Modes

# Forcing Shape
x0,κ  = ForcingParams()
ψ     = zeros(ComplexF64,Nby2)
SetForcingShape!(ψ,Bg,xg,x0,1.0,κ)
SEM1D.SEM_SetBC!(ψ,Inp.lbc,Inp.rbc)

ax2.plot(xg,real.(ψ) ,linewidth=2,linestyle="-", color=cm(2),label=L"\mathfrak{R}(ψ)")
ax2.plot(xg,imag.(ψ) ,linewidth=2,linestyle="--",color=cm(2),label=L"\mathfrak{Im}(ψ)")

vzro  = 0.0*ψ
Lθ    = zeros(ComplexF64,N,h)
f1    = 1.0/(sqrt(2.0))*[ψ;ψ]
f2    = [ψ;  vzro]
f3    = [vzro;ψ]
f4    = [ψ;  vzro]
f5    = [vzro;ψ]
f6    = [ψ; conj.(ψ)]
f7    = [conj.(ψ); ψ]


# Lθ    = [f1 f2 f3 f4 f5]
Lθ    = [f2 f3]
#Lθ    = [f6 f7]
#λh    = zeros(ComplexF64,h)
#λh    = [0.0im; im; -im; 2.3im; -2.3im]
#Lθ[:,1] = copy(f2)
if ifresonant
  λh    = [1.0im; -1.0im;]
else
  λh    = [2.3im; -2.3im;]
#  λh    = [0.7im; -0.7im;]
#  λh    = [2.3im]
end
Λh    = diagm(λh)

ΓH    = zeros(ComplexF64,n,h)
for i in 1:h
  for j in 1:n
    if abs(λc[j] - λh[i]) < 1.0e-12
      ΓH[j,i]   = W[:,j]'*(Bg2.*Lθ[:,i])
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
Ip    = zeros(ComplexF64,p,p)
Rp    = zeros(ComplexF64,N,p)
ind1  = 1:Nby2
ind2  = Nby2+1:N
for i in 1:p
  r  = Bg2M*copy(Lν[:,i])
  DQ = Matrix{ComplexF64}(I,Nby2,Nby2)
  CQ = Matrix{ComplexF64}(I,Nby2,Nby2)
  for j in 1:n
    if abs(λc[j] - λp[i]) < 1.0e-12
      println("Resonant λp: $i, $j")
      r  .= r .- V[:,j]*ΓP[j,i]
      DQ .= DQ - V[ind1,j]*(W[ind1,j].*Bg)'
      CQ .= CQ - V[ind2,j]*(W[ind2,j].*Bg)'
    end
  end  
  Rp[:,i]         = copy(r)
  
  @views SEM1D.SEM_SetBC!(r[ind1],Inp.lbc,Inp.rbc)
  @views SEM1D.SEM_SetBC!(r[ind2],Inp.lbc,Inp.rbc)
 
  ω               = λp[i]
  Res1            = DQ*(ω*I - OPg)*DQ
  vp1             = copy(r[ind1])
  gmres!(vp1,Res1,r[ind1])
  Vp[ind1,i]      = copy(vp1)

  
  Res2            = CQ*(ω*I - OPCg)*CQ
  vp2             = copy(r[ind2])
  gmres!(vp2,Res2,r[ind2])
  Vp[ind2,i]      = copy(vp2)

  vnorm = sqrt(abs(Vp[:,i]'*(Bg2.*Vp[:,i])))

  Ip[i,i] = ComplexF64(1)
 
  # Plot Mode
  if (emodeplot) && vnorm > 0.0
    j = n+i
    # ax2.plot(xg,real.(vp1),linewidth=2,linestyle="-", color=cm(j-1),label=L"\mathfrak{R}(ϕ_{%$j})")
    # ax2.plot(xg,imag.(vp1),linewidth=2,linestyle="--",color=cm(j-1),label=L"\mathfrak{Im}(ϕ_{%$j})")
  end

end  


# Forcing Modes
Vh    = zeros(ComplexF64,N,h)
Rh    = zeros(ComplexF64,N,h)
Ih    = zeros(ComplexF64,h,h)
# ind1  = 1:Nby2
# ind2  = Nby2+1:N
for i in 1:h
  r  = Bg2M*copy(Lθ[:,i])
  DQ = Matrix{ComplexF64}(I,Nby2,Nby2)
  CQ = Matrix{ComplexF64}(I,Nby2,Nby2)
  resonance = false
  for j in 1:n
    if abs(λc[j] - λh[i]) < 1.0e-12
      println("Resonant λh: $i, $j")
      resonance = true
      r  .= r .- V[:,j]*ΓH[j,i]
      DQ .= DQ - V[ind1,j]*(W[ind1,j].*Bg)'
      CQ .= CQ - V[ind2,j]*(W[ind2,j].*Bg)'
    end
  end  
  Rh[:,i]         = copy(r)
 
  @views SEM1D.SEM_SetBC!(r[ind1],Inp.lbc,Inp.rbc)
  @views SEM1D.SEM_SetBC!(r[ind2],Inp.lbc,Inp.rbc)
 
  ω               = λh[i]
  if (resonance)
    Res1          = DQ*(ω*BgM - OPg)*DQ
  else
    Res1          = (ω*BgM - OPg)
  end  
  vh1             = copy(r[ind1])
  @views gmres!(vh1,Res1,r[ind1])
  Vh[ind1,i]      = copy(vh1)

  if (resonance)
    Res2          = CQ*(ω*BgM - OPCg)*CQ
  else
    Res2          = (ω*BgM - OPCg)
  end  
  # Res2            = CQ*(ω*I - OPCg)*CQ
  vh2             = copy(r[ind2])
  @views gmres!(vh2,Res2,r[ind2])
  Vh[ind2,i]      = copy(vh2)

  vnorm = sqrt(abs(Vh[:,i]'*(Bg2.*Vh[:,i])))

  Ih[i,i] = ComplexF64(1)

  # Plot Mode
  if (emodeplot) && imag.(λh[i]) >= 0.0
    j = n+p+i
    ax2.plot(xg,real.(vh1),linewidth=2,linestyle="-", color=cm(j-1),label=L"\mathfrak{R}(ϕ_{%$j})")
    ax2.plot(xg,imag.(vh1),linewidth=2,linestyle="--",color=cm(j-1),label=L"\mathfrak{Im}(ϕ_{%$j})")
  end
end  

Vext  = [V  Vp  Vh]


ZeroPN = zeros(ComplexF64,p,n)
ZeroHN = zeros(ComplexF64,h,n)
ZeroPP = zeros(ComplexF64,p,p)
ZeroHH = zeros(ComplexF64,h,h)
ZeroPH = zeros(ComplexF64,p,h)
ZeroHP = zeros(ComplexF64,h,p)

Vhat  = [V       Vp     Vh;
         ZeroPN  Ip     ZeroPH;
         ZeroHN  ZeroHP Ih]

if (emodeplot)
  ax2.legend(ncols=4,fontsize=Grh.lgfs)
else  
  ax2.legend(ncols=3,fontsize=Grh.lgfs)
end  


# Extended Adjoint Tangent Space
#-------------------------------------------------- 

# Adjoint Parameter modes
Wp    = zeros(ComplexF64,N,p)
Ip    = zeros(ComplexF64,p,p)
Zp    = Lν'*W
for j in 1:p
  for i in 1:n
    if abs(λc[i]' - λp[j]') < 1.0e-12
      println("Resonant λp*: $i, $j")
      Zp[j,i] = 0.0
    else
      Zp[j,i] = Zp[j,i]/(λc[i]' - λp[j]')
    end
  end
  Ip[j,j] = ComplexF64(1)
end  

# Adjoint Forcing modes
Wh    = zeros(ComplexF64,N,h)
Ih    = zeros(ComplexF64,h,h)
Zh    = Lθ'*(Bg2M*W)
for j in 1:h
  for i in 1:n
    if abs(λc[i]' - λh[j]') < 1.0e-12
      println("Resonant λh*: $i, $j")
      Zh[j,i] = 0.0
    else
      Zh[j,i] = Zh[j,i]/(λc[i]' - λh[j]')
    end
  end
  Ih[j,j] = ComplexF64(1)
end  

Wext  = [W  Wp  Wh]

What  = [W   Wp     Wh;
         Zp  Ip     ZeroPH;
         Zh  ZeroHP Ih]

ZeroPN      = zeros(ComplexF64,p,ndof*2)
ZeroHN      = zeros(ComplexF64,h,ndof*2)
ZeroNN      = zeros(ComplexF64,ndof,ndof)

OPg2        = [OPg      ZeroNN;
               ZeroNN   OPCg]

Lhat        = [OPg2     Lν          Lθ;
               ZeroPN   Λp          ZeroPH;
               ZeroHN   ZeroHP      Λh]


Bghat = [Bg2; ones(ComplexF64,p); ones(ComplexF64,h)]

BiOrtho = What'*(diagm(Bghat)*Vhat)

println("Extended Tangent Space Done.")

println("Extended Tangent Space using Arnoldi.")
# Stepper-Arnoldi
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
vlen              = ndof+1
nev               = 2
ekryl             = 15  
lkryl             = nev + ekryl 
ngs               = 2
bsize             = 1
outer_iterations  = 100
tol               = 1.0e-12
EArnInp           = StepperArnoldi.ArnoldiInput(ifarnoldi,ifverbose,vlen,nev,ekryl,lkryl,ngs,bsize,outer_iterations,tol)

# Direct
EOPg                 = spzeros(ComplexF64,EArnInp.vlen,EArnInp.vlen)
EOPg[1:ndof,1:ndof]  = OPg
EOPg[1:ndof,ndof+1]  = Bg.*ψ
EOPg[ndof+1,ndof+1]  = λh[1]
EBg                  = [Bg[:]; 1.0]
EArnDir              = StepperArnoldi.StepArn( EOPg,EBg,EStpInp,EArnInp,Inp.lbc,Inp.rbc)

# Adjoint
EAOPg                = spzeros(ComplexF64,EArnInp.vlen,EArnInp.vlen)
EAOPg[1:ndof,1:ndof] = AOPg
EAOPg[ndof+1,1:ndof] = (Bg.*ψ)'
EAOPg[ndof+1,ndof+1] = λh[1]'
EArnAdj              = StepperArnoldi.StepArn( EAOPg,EBg,EStpInp,EArnInp,Inp.lbc,Inp.rbc)

ind7        = argmin(abs.(EArnDir.evals .- λh[1]))
θt          = EArnDir.evecs[ndof+1,ind7]
v7          = EArnDir.evecs[:,ind7]./θt 
ax2.plot(xg,real.(v7[ind1]),linewidth=2,linestyle="-", color=cm(7-1),label=L"\mathfrak{R}(ϕ_{7})")
ax2.plot(xg,imag.(v7[ind1]),linewidth=2,linestyle="--",color=cm(7-1),label=L"\mathfrak{Im}(ϕ_{7})")

ind8        = argmin(abs.(EArnAdj.evals .- λc[1]'))
θt          = EArnAdj.evecs[ndof+1,ind8]
θt2         = θt/Zh[1,1]         
v8          = EArnAdj.evecs[:,ind8]# ./θt2
#@views renormalize_evecs!(V[ind1,1],v8[ind1],Bg)
#ax2.plot(xg,real.(v8[ind1]),linewidth=1,linestyle="-", color=cm(8-1),label=L"\mathfrak{R}(χ_{7})")
#ax2.plot(xg,imag.(v8[ind1]),linewidth=1,linestyle="--",color=cm(8-1),label=L"\mathfrak{Im}(χ_{7})")

vtmp        = zeros(ComplexF64,ndof+1)
vtmp[ind1]  = copy(V[ind1])
wtmp        = zeros(ComplexF64,ndof+1)
wtmp        = copy(v8)
Bg3         = zeros(Float64,ndof+1)
Bg3[ind1]   = copy(Bg)
Bg3[ndof+1] = 1.0
@views renormalize_evecs!(vtmp,wtmp,Bg3)
ax2.plot(xg,real.(wtmp[ind1]),linewidth=1,linestyle="-", color=cm(8-1),label=L"\mathfrak{R}(χ_{7})")
ax2.plot(xg,imag.(wtmp[ind1]),linewidth=1,linestyle="--",color=cm(8-1),label=L"\mathfrak{Im}(χ_{7})")

Lθ1         = Lθ[ind1,1]
zzz1        = Lθ1'*wtmp[ind1]/(λc[1]' - λh[1]')
zzz2        = Lθ1'*(BgM*wtmp[ind1])/(λc[1]' - λh[1]')


if (emodeplot)
  ax2.legend(ncols=4,fontsize=Grh.lgfs)
else  
  ax2.legend(ncols=3,fontsize=Grh.lgfs)
end  

Vhat[ind1,5] = v7[ind1]

BiOrtho2 = What'*(diagm(Bghat)*Vhat)

println("Extended Tangent Space (Arnoldi) Done.")














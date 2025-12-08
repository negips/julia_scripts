# Testing the module

# Extending the tangent space
#---------------------------------------------------------------------- 

ifresonant = true

Nby2  = ArnInp.vlen
N     = Nby2*2
n     = 2
p     = 2
h     = 2
m     = n+p+h

Bg2   = [Bg; Bg]

# Parameter modes
Lν    = zeros(ComplexF64,N,p)
ΓP    = zeros(ComplexF64,n,p)
λp    = zeros(ComplexF64,p)
for i in 1:p
  for j in 1:n
    if abs(λc[j] - λp[i]) < 1.0e-12
      ΓH[j,i]   = W[:,j]'*(Bg2.*Lν[:,i])
    end
  end
end

# Forcing Modes

# Forcing Shape
x0    = ForcingLocation()
# ψ     = ForcingShape(Bg,xg,x0,1.0)
ψ     = zeros(ComplexF64,Nby2)
SetForcingShape!(ψ,Bg,xg,x0,1.0)

ax2.plot(xg,real.(ψ) ,linewidth=2,linestyle="-", color=cm(2),label=L"\mathfrak{R}(ψ)")
ax2.plot(xg,imag.(ψ) ,linewidth=2,linestyle="--",color=cm(2),label=L"\mathfrak{Im}(ψ)")

vzro  = 0.0*ψ
Lθ    = zeros(ComplexF64,N,h)
f1    = 1.0/(sqrt(2.0))*[ψ;ψ]
f2    = [ψ;  vzro]
f3    = [vzro;ψ]
f4    = [ψ;  vzro]
f5    = [vzro;ψ]
# Lθ    = [f1 f2 f3 f4 f5]
Lθ    = [f2 f3]
#λh    = zeros(ComplexF64,h)
#λh    = [0.0im; im; -im; 2.3im; -2.3im]
if ifresonant
  λh    = [1.0im; -1.0im;]
else
  λh    = [1.3im; -1.3im;]
#  λh    = [0.7im; -0.7im;]
end

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
Rp    = zeros(ComplexF64,N,p)
ind1  = 1:Nby2
ind2  = Nby2+1:N
for i in 1:p
  r  = copy(Lν[:,i])
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
end  


# Forcing Modes
Vh    = zeros(ComplexF64,N,h)
Rh    = zeros(ComplexF64,N,h)
ind1  = 1:Nby2
ind2  = Nby2+1:N
for i in 1:h
  r  = copy(Lθ[:,i])
  DQ = Matrix{ComplexF64}(I,Nby2,Nby2)
  CQ = Matrix{ComplexF64}(I,Nby2,Nby2)
  for j in 1:n
    if abs(λc[j] - λh[i]) < 1.0e-12
      println("Resonant λh: $i, $j")
      r  .= r .- V[:,j]*ΓH[j,i]
      DQ .= DQ - V[ind1,j]*(W[ind1,j].*Bg)'
      CQ .= CQ - V[ind2,j]*(W[ind2,j].*Bg)'
    end
  end  
  Rh[:,i]         = copy(r)
 
  @views SEM1D.SEM_SetBC!(r[ind1],Inp.lbc,Inp.rbc)
  @views SEM1D.SEM_SetBC!(r[ind2],Inp.lbc,Inp.rbc)
 
  ω               = λh[i]
  Res1            = DQ*(ω*I - OPg)*DQ
  vh1             = copy(r[ind1])
  @views gmres!(vh1,Res1,r[ind1])
  Vh[ind1,i]      = copy(vh1)

  Res2            = CQ*(ω*I - OPCg)*CQ
  vh2             = copy(r[ind2])
  gmres!(vh2,Res2,r[ind2])
  Vh[ind2,i]      = copy(vh2)

  vnorm = sqrt(abs(Vh[:,i]'*(Bg2.*Vh[:,i])))
end  

Vext  = [V  Vp  Vh]
ax2.legend(ncols=3)


println("Extended Tangent Space Done.")
















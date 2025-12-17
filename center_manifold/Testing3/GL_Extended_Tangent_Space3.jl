# Testing the module

# Extending the tangent space
#---------------------------------------------------------------------- 
println("Extended Tangent Space using Arnoldi.")

ifresonant = true
emodeplot  = true

Nby2  = ArnInp.vlen
N     = Nby2*2
n     = length(λc)
p     = 2
h     = 2
m     = n+p+h

Bg2   = [Bg; Bg]
Bg2M  = diagm(Bg2)

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
vlen              = ndof+1
nev               = 1
ekryl             = 15  
lkryl             = nev + ekryl 
eigshift          = 0.0 + 0.0im
ngs               = 2
bsize             = 1
outer_iterations  = 100
tol               = 1.0e-12
EArnInp           = StepperArnoldi.ArnoldiInput(ifarnoldi,ifverbose,ifeigshift,vlen,nev,ekryl,lkryl,eigshift,ngs,bsize,outer_iterations,tol)

# Extended Direct
EOPg                  = spzeros(ComplexF64,EArnInp.vlen,EArnInp.vlen)
EOPg[1:ndof,1:ndof]   = OPg
EOPg[1:ndof,ndof+1]   = zeros(ComplexF64,ndof)
EOPg[ndof+1,ndof+1]   = 0.0
EBg                   = [Bg[:]; 1.0]

# Extended Conjugated Direct
EOPCg                 = spzeros(ComplexF64,EArnInp.vlen,EArnInp.vlen)
EOPCg[1:ndof,1:ndof]  = OPCg
EOPCg[1:ndof,ndof+1]  = zeros(ComplexF64,ndof)
EOPCg[ndof+1,ndof+1]  = 0.0
EBCg                  = [Bg[:]; 1.0]
#-------------------------------------------------- 



# Parameter modes
Lν    = zeros(ComplexF64,N,p)
ΓP    = zeros(ComplexF64,n,p)
λp    = zeros(ComplexF64,p)
Λp    = diagm(λp)
for i in 1:p
  for j in 1:n
    if abs(λc[j] - λp[i]) < 1.0e-12
      ΓP[j,i]   = W[:,j]'*(Bg2.*Lν[:,i])
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

# for i in 1:h
#   f = copy(Lθ[:,i])
#   for ig in 1:ngs
#     α = zeros(T2,n)
#     for j in 1:n
#       if abs(λc[j] - λe) < 1.0e-12
#         α[j]  = W[:,j]'*(B.*f)
#         ΓH[j,i]  = ΓH[j,i] + W[:,j]'*(B.*f)
#       end
#     end
#     f .= f .- V*α
#   end
# end



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

  Vr1   = zeros(ComplexF64,Nby2+1,n)
  Vr2   = zeros(ComplexF64,Nby2+1,n)

  Wr1   = zeros(ComplexF64,Nby2+1,n)
  Wr2   = zeros(ComplexF64,Nby2+1,n)

  EOPg[1:ndof,1:ndof]  = OPg

  r  = copy(Lν[:,i])
  resonance = false
  nr = 0
  for j in 1:n
    if abs(λc[j] - λp[i]) < 1.0e-12
      println("Resonant λh: $i, $j")
      resonance = true
      nr = nr + 1
      r  .= r .- V[:,j]*ΓP[j,i]

      di = (nr-1)*Nby2+1 + 1    # destination index
      si = (j-1)*N + 1        # source index
      copyto!(Vr1,di,V,si,Nby2)
      copyto!(Wr1,di,W,si,Nby2)

      di = (nr-1)*Nby2+1 + 1    # destination index
      si = (j-1)*N + Nby2 + 1 # source index
      copyto!(Vr2,di,V,si,Nby2)
      copyto!(Wr2,di,W,si,Nby2)
     
    end
  end  
  @views SEM1D.SEM_SetBC!(r[ind1],Inp.lbc,Inp.rbc)
  @views SEM1D.SEM_SetBC!(r[ind2],Inp.lbc,Inp.rbc)
  Rp[:,i]               = copy(r)

  ω                     = λp[i]
  EArnInp.eigshift      = ω
  if norm(r[ind1]) > EArnInp.tol
    EOPg[1:ndof,ndof+1]   = Bg.*r[ind1]
    EOPg[ndof+1,ndof+1]   = ω

    if (resonance)
      Vr1_view    = view(Vr1,:,1:nr)
      Wr1_view    = view(Wr1,:,1:nr)
      EArnDir1    = StepperArnoldi.RestrictedStepArn( EOPg,EBg,Vr1_view,Wr1_view,EStpInp,EArnInp,Inp.lbc,Inp.rbc) 
    else
      EArnDir1    = StepperArnoldi.StepArn( EOPg,EBg,EStpInp,EArnInp,Inp.lbc,Inp.rbc) 
    end
    ii            = argmin(abs.(EArnDir1.evals .- ω))
    θt            = EArnDir1.evecs[ndof+1,ii]
    vp1           = EArnDir1.evecs[1:ndof,ii]./θt     # ensure extended variable == 1.0 
    Vp[ind1,i]    = copy(vp1)
  else
    Vp[ind1,i]    = zeros(ComplexF64,ndof)
  end

  # Conjugated
  if norm(r[ind2]) > EArnInp.tol
    EOPCg[1:ndof,ndof+1]  = Bg.*r[ind2]
    EOPCg[ndof+1,ndof+1]  = ω
    if (resonance)
      Vr2_view    = view(Vr2,:,1:nr)
      Wr2_view    = view(Wr2,:,1:nr)
      EArnDir2    = StepperArnoldi.RestrictedStepArn( EOPCg,EBCg,Vr2_view,Wr2_view,EStpInp,EArnInp,Inp.lbc,Inp.rbc) 
    else
      EArnDir2    = StepperArnoldi.StepArn( EOPCg,EBCg,EStpInp,EArnInp,Inp.lbc,Inp.rbc) 
    end
    ii            = argmin(abs.(EArnDir2.evals .- ω))
    θt            = EArnDir2.evecs[ndof+1,ii]
    vp2           = EArnDir2.evecs[1:ndof,ii]./θt     # ensure extended variable == 1.0 
    Vp[ind2,i]    = copy(vp2)
  else
    Vp[ind2,i]    = zeros(ComplexF64,ndof)
  end

  vnorm = sqrt(abs(Vp[:,i]'*(Bg2.*Vp[:,i])))

  Ip[i,i] = ComplexF64(1)

  # Plot Mode
  if (emodeplot) && vnorm > 0.0
    j = n+i
    ax2.plot(xg,real.(vp1),linewidth=2,linestyle="-", color=cm(j-1),label=L"\mathfrak{R}(ϕ_{%$j})")
    ax2.plot(xg,imag.(vp1),linewidth=2,linestyle="--",color=cm(j-1),label=L"\mathfrak{Im}(ϕ_{%$j})")
  end

end  


# Forcing Modes
Vh    = zeros(ComplexF64,N,h)
Rh    = zeros(ComplexF64,N,h)
Ih    = zeros(ComplexF64,h,h)
for i in 1:h

  Vr1   = zeros(ComplexF64,Nby2+1,n)
  Vr2   = zeros(ComplexF64,Nby2+1,n)

  Wr1   = zeros(ComplexF64,Nby2+1,n)
  Wr2   = zeros(ComplexF64,Nby2+1,n)

  EOPg[1:ndof,1:ndof]  = OPg

  r  = copy(Lθ[:,i])
  resonance = false
  nr        = 0
  for j in 1:n
    if abs(λc[j] - λh[i]) < 1.0e-12
      println("Resonant λh: $i, $j")
      resonance = true
      nr  = nr + 1
      r  .= r .- V[:,j]*ΓH[j,i]
 
      di = (nr-1)*Nby2+1 + 1    # destination index
      si = (j-1)*N + 1        # source index
      copyto!(Vr1,di,V,si,Nby2)
      copyto!(Wr1,di,W,si,Nby2)

      di = (nr-1)*Nby2+1 + 1    # destination index
      si = (j-1)*N + Nby2 + 1 # source index
      copyto!(Vr2,di,V,si,Nby2)
      copyto!(Wr2,di,W,si,Nby2)
    end
  end  
  @views SEM1D.SEM_SetBC!(r[ind1],Inp.lbc,Inp.rbc)
  @views SEM1D.SEM_SetBC!(r[ind2],Inp.lbc,Inp.rbc)
  Rh[:,i]               = copy(r)

  ω                     = λh[i]
  EArnInp.eigshift      = ω
  if norm(r[ind1]) > EArnInp.tol
    EOPg[1:ndof,ndof+1]   = Bg.*r[ind1]
    EOPg[ndof+1,ndof+1]   = ω

    if (resonance)
      Vr1_view    = view(Vr1,:,1:nr)
      Wr1_view    = view(Wr1,:,1:nr)
      EArnDir1    = StepperArnoldi.RestrictedStepArn( EOPg,EBg,Vr1_view,Wr1_view,EStpInp,EArnInp,Inp.lbc,Inp.rbc) 
    else
      EArnDir1    = StepperArnoldi.StepArn( EOPg,EBg,EStpInp,EArnInp,Inp.lbc,Inp.rbc) 
    end
    ii            = argmin(abs.(EArnDir1.evals .- ω))
    θt            = EArnDir1.evecs[ndof+1,ii]
    vh1           = EArnDir1.evecs[1:ndof,ii]./θt     # ensure extended variable == 1.0 
    Vh[ind1,i]    = copy(vh1)
  else
    Vh[ind1,i]    = zeros(ComplexF64,ndof)
  end

  # Conjugated
  if norm(r[ind2]) > EArnInp.tol
    EOPCg[1:ndof,ndof+1]  = Bg.*r[ind2]
    EOPCg[ndof+1,ndof+1]  = ω
    if (resonance)
      Vr2_view    = view(Vr2,:,1:nr)
      Wr2_view    = view(Wr2,:,1:nr)
      EArnDir2    = StepperArnoldi.RestrictedStepArn( EOPCg,EBCg,Vr2_view,Wr2_view,EStpInp,EArnInp,Inp.lbc,Inp.rbc) 
    else
      EArnDir2    = StepperArnoldi.StepArn( EOPCg,EBCg,EStpInp,EArnInp,Inp.lbc,Inp.rbc) 
    end
    ii            = argmin(abs.(EArnDir2.evals .- ω))
    θt            = EArnDir2.evecs[ndof+1,ii]
    vh2           = EArnDir2.evecs[1:ndof,ii]./θt     # ensure extended variable == 1.0 
    Vh[ind2,i]    = copy(vh2)
  else
    Vh[ind2,i]    = zeros(ComplexF64,ndof)
  end

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
#Zh    = Lθ'*(Bg2M*W)
Zh    = -Vh'*(Bg2M*W)
for j in 1:h
  for i in 1:n
    if abs(λc[i]' - λh[j]') < 1.0e-12
      println("Resonant λh*: $i, $j")
      Zh[j,i] = -Vh[:,j]'*(Bg2.*W[:,i])
    else
      # Zh[j,i] = Zh[j,i]/(λc[i]' - λh[j]')
      Zh[j,i] = -Vh[:,j]'*(Bg2.*W[:,i])
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

Lhat        = [OPg2     Bg2M*Lν   Bg2M*Lθ;
               ZeroPN   Λp        ZeroPH;
               ZeroHN   ZeroHP    Λh]

Bghat       = [Bg2; ones(ComplexF64,p); ones(ComplexF64,h)]

BiOrtho     = What'*(diagm(Bghat)*Vhat)

if (emodeplot)
  ax2.legend(ncols=4,fontsize=Grh.lgfs)
else  
  ax2.legend(ncols=3,fontsize=Grh.lgfs)
end  

println("Extended Tangent Space (Arnoldi) Done.")














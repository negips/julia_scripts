# Testing the module

# Extending the tangent space
#---------------------------------------------------------------------- 
println("Extended Tangent Space using Arnoldi.")
#---------------------------------------------------------------------- 
include("GL_Extended_Tangent_Space5.jl")
#----------------------------------------------------------------------

ifresonant = true
emodeplot  = true
restricted = false 

Nby2  = size(OPg,2)
N     = Nby2*2
n     = length(λc)
p     = 0
s     = npert
h     = 2
m     = nsys+npert+p+h

# Parameter Perturbation
Lν    = zeros(ComplexF64,N,p)
λν    = zeros(ComplexF64,p)

# System Perturbation
Lσ    = zeros(ComplexF64,N,s)
λσ    = zeros(ComplexF64,s)

# Forcing Modes
ξ,Lθ,λh     = GetExternalForcing(xg,Bg,ifresonant,Inp.lbc,Inp.rbc)
Λh          = diagm(λh)
j_cm        = Int(nsys)
ax2.plot(xg,real.(ξ) ,linewidth=2,linestyle="-", color=cm(j_cm),label=L"\mathfrak{R}(ξ)")
ax2.plot(xg,imag.(ξ) ,linewidth=2,linestyle="--",color=cm(j_cm),label=L"\mathfrak{Im}(ξ)")

λext        = zeros(ComplexF64,s+p+h)
#λext = [λσ; λν; λh]
Lext        = zeros(ComplexF64,N,s+p+h)
for i in 1:s
  for j in 1:N
    Lext[j,i] = Lσ[j,i]
  end
  λext[i] = λσ[i]
end
for i in 1:p
  for j in 1:N
    Lext[j,i+s] = Lν[j,i]
  end
  λext[s+i] = λν[i]
end
for i in 1:h
  for j in 1:N
    Lext[j,i+p+s] = Lθ[j,i]
  end
  λext[s+p+i] = λh[i]
end

if (ifmodepert)
  EM = GLExtendPertTangentSpace(OPg,OPCg,Bg,λSys,σ,VSys,WSys,λext,Lext,restricted,Inp.lbc,Inp.rbc)
else  
  EM = GLExtendTangentSpace(OPg,OPCg,Bg,λSys,VSys,WSys,λext,Lext,restricted,Inp.lbc,Inp.rbc)
end

for i in 1:length(λext)
  vnorm = norm(EM.Ve[ind1,i])
  # Plot Mode
  if (emodeplot) && vnorm > 0.0
    j = nsys + 1 + i
    ax2.plot(xg,real.(EM.Ve[ind1,i]),linewidth=2,linestyle="-", color=cm(j-1),label=L"\mathfrak{R}(ϕ_{%$j})")
    ax2.plot(xg,imag.(EM.Ve[ind1,i]),linewidth=2,linestyle="--",color=cm(j-1),label=L"\mathfrak{Im}(ϕ_{%$j})")
  end

  # # Conjugate Mode
  # vnorm2 = norm(EM.Ve[ind2,i])
  # # Plot Mode
  # if (emodeplot) && vnorm2 > 0.0
  #   j = Int(nsys/2) + 1 + i + 1
  #   ax2.plot(xg,real.(EM.Ve[ind2,i]),linewidth=2,linestyle="-", color=cm(j-1),label=L"\mathfrak{R}(ϕ_{%$j})")
  #   ax2.plot(xg,imag.(EM.Ve[ind2,i]),linewidth=2,linestyle="--",color=cm(j-1),label=L"\mathfrak{Im}(ϕ_{%$j})")
  # end

end  

if (ifmodepert)
  Vext        = [VSys EM.Ve]
  Wext        = [WSys EM.We]
  Γe          = EM.Γ
  Ze          = EM.Z
  Λe          = diagm(EM.λe)
  ΛSys        = diagm(λSys)
  Zero_ne_n   = zeros(ComplexF64,s+p+h,nsys)
else
  Vext        = [VSys EM.Ve]
  Wext        = [WSys EM.We]
  Γe          = EM.Γ
  Ze          = EM.Z
  Λe          = diagm(EM.λe)
  ΛSys        = diagm(λSys)
  Zero_ne_n   = zeros(ComplexF64,s+p+h,nsys)
end  


Khat  = [ΛSys      Γe;
         Zero_ne_n Λe]

Vhat  = [Vext;
         Zero_ne_n I]
What  = [Wext;
         Ze        I]
Bhat  = [Bg2; ones(eltype(Bg2),s+p+h)]

EBiOrtho = What'*diagm(Bhat)*Vhat

if (emodeplot)
  if h>0
    ax2.legend(ncols=4,fontsize=Grh.lgfs)
  else
    ax2.legend(ncols=3,fontsize=Grh.lgfs)
  end
end  
h2fname = "extended_eigenvectors.eps"
save_figure(h2,h2fname,figsave)  

println("Extended Tangent Space (Arnoldi) Done.")





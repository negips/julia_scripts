# Testing the module

# Extending the tangent space
#---------------------------------------------------------------------- 
println("Extended Tangent Space using Arnoldi.")
#---------------------------------------------------------------------- 
include("GL_Extended_Tangent_Space5.jl")
#----------------------------------------------------------------------

ifresonant = false
emodeplot  = true
restricted = false 

Nby2  = size(OPg,2)
N     = Nby2*2
n     = length(λc)
p     = 2
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
ψ,Lθ,λh     = GetExternalForcing(xg,Bg,Inp.lbc,Inp.rbc)
Λh          = diagm(λh)
ax2.plot(xg,real.(ψ) ,linewidth=2,linestyle="-", color=cm(nsys),label=L"\mathfrak{R}(ψ)")
ax2.plot(xg,imag.(ψ) ,linewidth=2,linestyle="--",color=cm(nsys),label=L"\mathfrak{Im}(ψ)")


# PE,HE = GLExtendTangentSpace2(OPg,OPCg,Bg,λc,V,W,λν,Lν,λh,Lθ,restricted,Inp.lbc,Inp.rbc)

λext = [λσ; λν; λh]
Lext = zeros(ComplexF64,N,p+s+h)
for i in 1:s
  for j in 1:N
    Lext[j,i] = Lσ[j,i]
  end
end
for i in 1:p
  for j in 1:N
    Lext[j,i+s] = Lν[j,i]
  end
end
for i in 1:h
  for j in 1:N
    Lext[j,i+p+s] = Lθ[j,i]
  end
end

if (ifmodepert)
  EM = GLExtendPertTangentSpace(OPg,OPCg,Bg,λSys,σ,VSys,WSys,λext,Lext,restricted,Inp.lbc,Inp.rbc)
else  
  EM = GLExtendTangentSpace(OPg,OPCg,Bg,λSys,VSys,WSys,λext,Lext,restricted,Inp.lbc,Inp.rbc)
end

# if (ifmodepert)
#   EM = GLExtendTangentSpace2OP(NewOPg,NewOPCg,Bg,λSys,VSys,WSys,λext,Lext,restricted,Inp.lbc,Inp.rbc)
# else
#   EM = GLExtendTangentSpace2(NewOPg,NewOPCg,Bg,λc,V,W,λext,Lext,restricted,Inp.lbc,Inp.rbc)
# end
# 
for i in 1:length(λext)
  vnorm = norm(EM.Ve[ind1,i])
  # Plot Mode
  if (emodeplot) && vnorm > 0.0
    j = nsys+i
    ax2.plot(xg,real.(EM.Ve[ind1,i]),linewidth=2,linestyle="-", color=cm(j-1),label=L"\mathfrak{R}(ϕ_{%$j})")
    ax2.plot(xg,imag.(EM.Ve[ind1,i]),linewidth=2,linestyle="--",color=cm(j-1),label=L"\mathfrak{Im}(ϕ_{%$j})")
  end
end  

if (ifmodepert)
  Vext        = [VSys EM.Ve]
  Wext        = [WSys EM.We]
  Γe          = EM.Γ
  Ze          = EM.Z
  Λe          = diagm(EM.λe)
  ΛSys        = diagm(λSys)
  Zero_ne_n   = zeros(ComplexF64,p+s+h,nsys)
else
  Vext        = [VSys EM.Ve]
  Wext        = [WSys EM.We]
  Γe          = EM.Γ
  Ze          = EM.Z
  Λe          = diagm(EM.λe)
  ΛSys        = diagm(λSys)
  Zero_ne_n   = zeros(ComplexF64,p+s+h,nsys)
end  


Khat  = [ΛSys      Γe;
         Zero_ne_n Λe]

Vhat  = [Vext;
         Zero_ne_n I]
What  = [Wext;
         Ze        I]
Bhat  = [Bg2; ones(eltype(Bg2),p+s+h)]

EBiOrtho = What'*diagm(Bhat)*Vhat

if (emodeplot)
  ax2.legend(ncols=4,fontsize=Grh.lgfs)
else  
  ax2.legend(ncols=3,fontsize=Grh.lgfs)
end  

println("Extended Tangent Space (Arnoldi) Done.")














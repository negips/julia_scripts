# Testing the module

#include("Module_SEM1D/SEM1D.jl")
#using .SEM1D

include("Module_StepperArnoldi/StepperArnoldi.jl")
#using .StepperArnoldi

#include("Module_CenterManifold/CenterManifold.jl")
#using .CenterManifold

#using LinearAlgebra
#using SparseArrays
#using Printf
#using PolynomialBases
#using IterativeSolvers
#using PyPlot

#---------------------------------------------------------------------- 


# Stepper-Arnoldi
#-------------------------------------------------- 
ifadjoint         = false
ifoptimal         = false
ifverbose         = false
verbosestep       = 2000
nsteps            = 2000
dt                = 2.5e-5
StpInp            = StepperArnoldi.StepperInput(ifadjoint,ifoptimal,ifverbose,verbosestep,nsteps,dt)

ifarnoldi         = true 
ifverbose         = false
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
cm    = get_cmap("tab10")

# Build Center-Manifold matrices
v1    = copy(ArnDir.evecs[:,1])
w1    = copy(ArnAdj.evecs[:,1])

fx0   = 6.0
j0    = argmin(abs.(xg .- fx0))
renormalize_evec!(v1,j0)
renormalize_evecs!(v1,w1,Bg)
v2    = conj.(v1)
w2    = conj.(w1)
λc    = [im; -im;]

# Forcing
# ψ     = exp.(-(xg .- fx0).^2) 
# ψn    = sqrt(ψ'*(Bg.*ψ))
# ψ    .= ψ./ψn

h2    = figure(num=2,figsize=[12.,9.])
ax2   = gca()
for i in 1:ArnInp.nev
  vtmp = ArnDir.evecs[:,i]
  wtmp = ArnAdj.evecs[:,i]
  renormalize_evec!(vtmp,j0)
  renormalize_evecs!(vtmp,wtmp,Bg)

  ax2.plot(xg,real.(vtmp),linewidth=2,linestyle="-", color=cm(i-1),label=L"\mathfrak{R}(ϕ_{%$i})")
  ax2.plot(xg,imag.(vtmp),linewidth=2,linestyle="--",color=cm(i-1),label=L"\mathfrak{Im}(ϕ_{%$i})")

  ax2.plot(xg,real.(wtmp),linewidth=1,linestyle="-", color=cm(i+nev-1),label=L"\mathfrak{R}(χ_{%$i})")
  ax2.plot(xg,imag.(wtmp),linewidth=1,linestyle="--",color=cm(i+nev-1),label=L"\mathfrak{Im}(χ_{%$i})")
end

# ax2.plot(xg,real.(v1),linewidth=2,linestyle="-", color=cm(0),label=L"\mathfrak{R}(ϕ)")
# ax2.plot(xg,imag.(v1),linewidth=2,linestyle="--",color=cm(0),label=L"\mathfrak{Im}(ϕ)")
# ax2.plot(xg,real.(w1),linewidth=2,linestyle="-", color=cm(1),label=L"\mathfrak{R}(χ)")
# ax2.plot(xg,imag.(w1),linewidth=2,linestyle="--",color=cm(1),label=L"\mathfrak{Im}(χ)")
#ax2.plot(xg,real.(ψ) ,linewidth=2,linestyle="-", color=cm(2),label=L"\mathfrak{R}(ψ)")
#ax2.plot(xg,imag.(ψ) ,linewidth=2,linestyle="--",color=cm(2),label=L"\mathfrak{Im}(ψ)")

vzro  = 0.0*v1

V     = [v1       vzro;
         vzro     v2]
W     = [w1       vzro;
         vzro     w2]

println("Standard Tangent Space Done.")
















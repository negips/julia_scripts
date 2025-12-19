# Standard Tangent Space

#include("Module_SEM1D/SEM1D.jl")
#using .SEM1D

include("../Module_StepperArnoldi/StepperArnoldi.jl")
#using .StepperArnoldi

#---------------------------------------------------------------------- 
screen = 1
Grh    = setgraphics(screen)

# Stepper-Arnoldi
#-------------------------------------------------- 
ifadjoint         = false
ifoptimal         = false
ifverbose         = false
verbosestep       = 500
nsteps            = 500
dt                = 1.0e-4
StpInp            = StepperArnoldi.StepperInput(ifadjoint,ifoptimal,ifverbose,verbosestep,nsteps,dt)
#StpInp      = Set_StepperParams()

ifarnoldi         = true 
ifverbose         = false
ifeigshift        = false
vlen              = ndof
nev               = 2
ekryl             = 15  
lkryl             = nev + ekryl 
eigshift          = 0.0 + 1.0im
ngs               = 2
bsize             = 1
outer_iterations  = 100
tol               = 1.0e-12
ArnInp            = StepperArnoldi.ArnoldiInput(ifarnoldi,ifverbose,ifeigshift,vlen,nev,ekryl,lkryl,eigshift,ngs,bsize,outer_iterations,tol)
#ArnInp      = Set_ArnoldiParams()

Ω     = SEM1D.GLAnalyticalSpectra(δ)
if (StpInp.ifadjoint)
  Ω   = conj.(Ω)
end  

#figz1   = [10.,7.]

close("all")
h1    = figure(num=1,figsize=Grh.figsz1);
ax1   = gca()
pl1   = ax1.plot(imag.(Ω),real.(Ω),linestyle="none",marker="o",markersize=1.5*Grh.mksz,markerfacecolor="none",markeredgewidth=3)
ax1.set_ylim([-2.75,0.5])
if (StpInp.ifadjoint)
  ax1.set_xlim([-1.80,0.0])
else  
  ax1.set_xlim([0.0,1.80])
end  

ArnDir      = StepperArnoldi.StepArn( OPg,Bg,StpInp,ArnInp,Inp.lbc,Inp.rbc)
ArnAdj      = StepperArnoldi.StepArn(AOPg,Bg,StpInp,ArnInp,Inp.lbc,Inp.rbc)

pl3   = ax1.plot(imag.(ArnDir.evals),real.(ArnDir.evals),linestyle="none",marker="o",markersize=Grh.mksz)
cm    = get_cmap("tab10")

# Build Center-Manifold matrices
id    = argmin(abs.(ArnDir.evals .- Ω[1]))
v1    = copy(ArnDir.evecs[:,id])
ia    = argmin(abs.(ArnAdj.evals .- Ω[1]'))
w1    = copy(ArnAdj.evecs[:,ia])

fx0   = 6.0
j0    = argmin(abs.(xg .- fx0))
renormalize_evec!(v1,j0)
renormalize_evecs!(v1,w1,Bg)
v2    = conj.(v1)
w2    = conj.(w1)
λc    = [1.0im; -1.0im]
Λc    = diagm(λc)
n     = length(λc)

# Plot (normalized) Eigenvectors
h2    = figure(num=2,figsize=Grh.figsz2)
ax2   = gca()
for i in 1:ArnInp.nev
  vtmp = ArnDir.evecs[:,i]
  wtmp = ArnAdj.evecs[:,i]
  renormalize_evec!(vtmp,j0)
  renormalize_evecs!(vtmp,wtmp,Bg)

  ax2.plot(xg,real.(vtmp),linewidth=2,linestyle="-", color=cm(i-1),label=L"\mathfrak{R}(ϕ_{%$i})")
  ax2.plot(xg,imag.(vtmp),linewidth=2,linestyle="--",color=cm(i-1),label=L"\mathfrak{Im}(ϕ_{%$i})")

  ax2.plot(xg,real.(wtmp),linewidth=1,linestyle="-", color=cm(i+ArnInp.nev-1),label=L"\mathfrak{R}(χ_{%$i})")
  ax2.plot(xg,imag.(wtmp),linewidth=1,linestyle="--",color=cm(i+ArnInp.nev-1),label=L"\mathfrak{Im}(χ_{%$i})")
end
ax2.set_xlabel(L"x",fontsize=Grh.lafs)
ax2.set_ylabel(L"A",fontsize=Grh.lafs)


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

Nby2  = size(OPg,2)
N     = 2*Nby2
ind1  = 1:Nby2
ind2  = Nby2+1:N
Bg2   = [Bg; Bg]
Bg2M  = diagm(Bg2)

println("Standard Tangent Space Done.")
















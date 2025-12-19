
include("../Module_StepperArnoldi/StepperArnoldi.jl")

# Stepper-Arnoldi
#-------------------------------------------------- 
ifadjoint         = false
ifoptimal         = false
ifverbose         = false
verbosestep       = 500
nsteps            = 500
dt                = 1.0e-4
StpInp2           = StepperArnoldi.StepperInput(ifadjoint,ifoptimal,ifverbose,verbosestep,nsteps,dt)
#StpInp      = Set_StepperParams()

ifarnoldi         = true 
ifverbose         = false
ifeigshift        = true
vlen              = ndof+1
nev               = 1
ekryl             = 15  
lkryl             = nev + ekryl 
eigshift          = λext[5]
ngs               = 2
bsize             = 1
outer_iterations  = 100
tol               = 1.0e-12
ArnInp2           = StepperArnoldi.ArnoldiInput(ifarnoldi,ifverbose,ifeigshift,vlen,nev,ekryl,lkryl,eigshift,ngs,bsize,outer_iterations,tol)

#
λf                = λext[5]
ff                = copy(Lext[ind1,5])
ifres             = fill(false,nsys)

DArn = StepperArnoldi.REPStepArn(OPg,Bg,σ,VSys[ind1,:],WSys[ind1,:],ifres,ff,λf,StpInp2,ArnInp2,Inp.lbc,Inp.rbc);

bb1 = [Bg; 1.0]
pl1 = ax2.plot(xg,real.(DArn.evecs[ind1]./DArn.evecs[ndof+1]),linewidth=2,linestyle="-", color="black",label=L"\mathfrak{R}(ϕ_{20})")
pl2 = ax2.plot(xg,imag.(DArn.evecs[ind1]./DArn.evecs[ndof+1]),linewidth=2,linestyle="--",color="black",label=L"\mathfrak{Im}(ϕ_{20})")

# Test REPLx
xt1         = copy(DArn.evecs[:,1])
dxt1        = StepperArnoldi.REPLx(xt1,OPg,Bg,VSys[ind1,:],WSys[ind2,:],σ,fill(false,nsys),ff,λf)
dxt1_B      = 0.0*dxt1
for i in 1:ndof
  dxt1_B[i] = dxt1[i]/Bg[i]
end
dxt1_B[ndof+1] = dxt1[ndof+1]/1.0

dxt1_B_xt1     = dxt1_B./xt1

println("Done.")

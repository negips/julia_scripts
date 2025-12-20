
include("../Module_StepperArnoldi/StepperArnoldi.jl")

# Stepper-Arnoldi
#-------------------------------------------------- 
ifadjoint         = false
ifoptimal         = false
ifverbose         = false
verbosestep       = 1000
nsteps            = 1000
dt                = 5.0e-5
StpInp2           = StepperArnoldi.StepperInput(ifadjoint,ifoptimal,ifverbose,verbosestep,nsteps,dt)
#StpInp      = Set_StepperParams()

ifarnoldi         = true 
ifverbose         = false
ifeigshift        = true
vlen              = ndof+1
nev               = 3
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

jj    = argmin(abs.(DArn.evals .- λf))
θr    = DArn.evecs[ndof+1,jj]

pl1 = ax2.plot(xg,real.(DArn.evecs[ind1,jj]./θr),linewidth=2,linestyle="-", color="black",label=L"\mathfrak{R}(ϕ_{20})")
pl2 = ax2.plot(xg,imag.(DArn.evecs[ind1,jj]./θr),linewidth=2,linestyle="--",color="black",label=L"\mathfrak{Im}(ϕ_{20})")


# Check
chking = (λf*diagm(Bg) .- OPg)*DArn.evecs[ind1,jj]/θr .- Bg.*ff

sol    = zeros(ComplexF64,ndof)
Bff    = Bg.*ff
ResOP  = (λf*diagm(Bg) .- OPg)

gmres!(sol,ResOP,Bff,abstol=1.0e-12)
pl3 = ax2.plot(xg,real.(sol),linewidth=2,linestyle="-", color="green",label=L"\mathfrak{R}(ϕ_{21})")
pl4 = ax2.plot(xg,imag.(sol),linewidth=2,linestyle="--",color="green",label=L"\mathfrak{Im}(ϕ_{21})")
# Check
chking2 = (λf*diagm(Bg) .- OPg)*sol .- Bff

Lbig    = zeros(ComplexF64,ndof+1,ndof+1)
Lbig[1:ndof,1:ndof] = copy(Matrix(OPg))
Lbig[1:ndof,ndof+1] = Bg.*ff
Lbig[ndof+1,ndof+1] = λf

Feig  = eigen(Lbig,diagm(bb1))
jj2   = argmin(abs.(Feig.values .- λf))
vv3   = Feig.vectors[:,jj2]
θr2   = vv3[ndof+1]
vv4   = vv3[1:ndof]./vv3[end]
pl5   = ax2.plot(xg,real.(vv4),linewidth=1,linestyle="--", color="yellow",label=L"\mathfrak{R}(ϕ_{22})")
pl6   = ax2.plot(xg,imag.(vv4),linewidth=1,linestyle="-.",color="yellow",label=L"\mathfrak{Im}(ϕ_{22})")



println("Done.")




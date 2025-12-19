
include("../Module_StepperArnoldi/StepperArnoldi.jl")

ArnInp.nev        = 1
ArnInp.vlen       = ndof + 1
ArnInp.ifeigshift = true
ArnInp.eigshift   = λext[5]
λf                = λext[5]
ff                = copy(Lext[ind1,5])
ifres             = fill(false,nsys)

DArn = StepperArnoldi.REPStepArn(OPg,Bg,σ,VSys[ind1,:],WSys[ind1,:],ifres,ff,λf,StpInp,ArnInp,Inp.lbc,Inp.rbc);

bb1 = [Bg; 1.0]
pl1 = ax2.plot(xg,real.(DArn.evecs[ind1]),linewidth=2,linestyle="-", color="black",label=L"\mathfrak{R}(ϕ_{20})")
pl2 = ax2.plot(xg,imag.(DArn.evecs[ind1]),linewidth=2,linestyle="--",color="black",label=L"\mathfrak{Im}(ϕ_{20})")

# Test REPLx


println("Done.")

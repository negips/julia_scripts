println("Newton interface")

using PolynomialBases
using PyPlot,PyCall
using LinearAlgebra
# using IterativeSolvers

include("$SRC/Dealias.jl")

include("custom_params.jl")
include("MyGMRES.jl")
include("newton_init.jl")
include("newton_op.jl")

X = Geom.xm1[:];
X[end] = X[1]
QTX = QT*(X.*vimult)

# No dynamic phase here
ifdynplot         = false
ifplot            = iffldplot || ifphplot

C                 = 1.0       # Front Velocity
Rhs               = 0.0*fld

dotfy(y)          = F(0.0,y)
dotfx(x)          = F(x,0.0)
Gradyf(y)         = GradFX(y,pars.fcy)/ϵ
Gradxf(x)         = GradFX(x,pars.fcx)/ϵ


Lkryl             = 20
VKryl             = zeros(VT,ndof,Lkryl+1)
tol               = 1.0e-8
maxoit            = 200
sol               = 0.0*fld[:,2]
δ                 = 1.0E-03

mask              = fill(1.0,ndof)

ub                = find_zero(dotfy,9.0) 
fld[1,2]          = ub

for i in 1:nsteps
  global fld,sol,C,mask
  global VKryl
  global pl

  # Boundary conditions: 
  # 0.0 at the right Boundary
  # Upper branch solution at the left boundary
  #fld[end,2]    = 0.0
  yguess        = 8.0
  upper_branch  = find_zero(dotfy,yguess)
  #fld[1,2]      = upper_branch

  # Masks are somewhat irrelevant for spectral approximations

  # Initial Residual
  j     = 2
  resid = newton_resid_spectral(fld[:,j],Basis,C,dotfy,Lg,Cg,M2N,N2M,upper_branch,0.0)
  
  # Boundary values are Dirichlet
  fnodal    = M2N*fld[:,j]
  grad_diag = Gradyf(fnodal)
  
  opg(x)    = newton_Lop_spectral(x,Basis,C,grad_diag,Lg,Cg,M2N,N2M)

  MyGMRES!(sol,resid,opg,VKryl,tol,maxoit)

  fld[:,j] .+= δ*sol

  resid = newton_resid_spectral(fld[:,j],Basis,C,dotfy,Lg,Cg,M2N,N2M,upper_branch,0.0)

  res2 = norm(resid)

  if mod(i,verbosestep)==0
    println("$i Residual: $res2")
  end

  if (res2 < tol)
    local pll = plot(x2,intpd*fld[:,j],color="black",linewidth=2);
    break
  end

  if (mod(i,plotupd)==0)
    for k in 1:nflds
      if (plotfldi[k])
        pl[k][1].remove()
        pl[k] = ax2.plot(x2,intpd*fld[:,k],color=cm(k-1));
        h2.show()
      end
    end  
  end  

end

if !@isdefined(h3)
  h3 = figure(num=3)
  pl3 = plot(x2,intpd*fld[:,2],color=cm(4),linewidth=2);
else
 figure(h3)
 h3.clf()
 h3.show()
 pl3 = plot(x2,intpd*fld[:,2],color=cm(4),linewidth=2);
end 





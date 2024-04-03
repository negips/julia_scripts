println("Newton interface")

using PolynomialBases
using PyPlot,PyCall
using LinearAlgebra
using IterativeSolvers

include("$SRC/Dealias.jl")

include("MyGMRES.jl")
include("newton_init.jl")
include("newton_op.jl")

X = Geom.xm1[:];
X[end] = X[1]
QTX = QT*(X.*vimult)

# No dynamic phase here
ifdynplot         = false
ifplot            = iffldplot || ifphplot

C                 = 4.0       # Front Velocity
Rhs               = 0.0*fld

dotfy(y)          = F(0.0,y)
dotfx(x)          = F(x,0.0)
Gradyf(y)         = GradFX(y,pars.fcy)/ϵ
Gradxf(x)         = GradFX(x,pars.fcx)/ϵ


Lkryl             = 40
VKryl             = zeros(VT,ndof,Lkryl+1)
tol               = 1.0e-8
maxoit            = 200
sol               = 0.0*fld[:,2]
δ                 = 0.5

mask              = fill(1.0,ndof)

for i in 1:nsteps
  global fld,sol,C,mask
  global VKryl

  # dotfld    = Flow(fld[:,1],fld[:,2])

  # Boundary conditions: 
  # 0.0 at the right Boundary
  # Upper branch solution at the left boundary
  fld[end,2]    = 0.0
  yguess        = 8.0
  upper_branch  = find_zero(dotfy,yguess)
  fld[1,2]      = upper_branch

  mask[1]       = 0.0
  mask[end]     = 0.0

  # Initial Residual
  j          = 2
  dotfld     = dotfy(fld[:,j]);
  lapfld     = Lg*fld[:,j];
  convfld    = Cg*fld[:,j];

  # global resid   = -mask.*(lapfld .+ C*convfld .+ Bg.*dotfld);
  resid = newton_resid(fld[:,j],C,dotfy,Lg,Cg,Bg,0.0,0.0)
  
  # Boundary values are Dirichlet
  grad_diag         = Gradyf(fld[:,2])
  
  # global opg(x)     = mask.*(Lg*x .+ C*Cg*x .+ Bg.*grad_diag.*x)
  opg(x)     = newton_Lop(x,C,grad_diag,Lg,Cg,Bg)

  MyGMRES!(sol,resid,opg,VKryl,tol,maxoit)


  fld[:,j] .+= δ*sol

  # dotfld     = dotfy(fld[:,j]);
  # lapfld     = Lg*fld[:,j];
  # convfld    = Cg*fld[:,j];
  # resid      = -mask.*(lapfld .+ C*convfld .+ Bg.*dotfld);

  resid = newton_resid(fld[:,j],C,dotfy,Lg,Cg,Bg,0.0,0.0)

  res2 = norm(resid)

  println("$i Residual: $res2")
  if (res2 < tol)
    local pl = plot(Geom.xm1[:],Q*fld[:,2],color="black",linewidth=2);
    break
  end

  local pl = plot(Geom.xm1[:],Q*fld[:,2],color=cm(2));

end

# close("all")

#ax2.clear()
#h2.clf()
#l=2
#pl[1] = plot(Geom.xm1[:],Q*fld[:,2],color=cm(l-1));

# sol2 = copy(sol)
# gmres!(sol2,M,resid)
# ax2.plot(Geom.xm1[:],Q*sol2,color=cm(l));

h3 = figure(num=3)
#figure(h3)
pl3 = plot(Geom.xm1[:],Q*fld[:,2],color=cm(3),linewidth=2);





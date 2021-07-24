println("Main interface for 1D SEM")

using PolynomialBases
using LinearAlgebra
# using UnicodePlots
#using Plots

# Include the function files
include("sem_geom.jl")

include("AssembleMatrix.jl")

include("AssembleMatrixLesshafft.jl")

include("sem_init_ref.jl")

include("Sem_QQT.jl")

Geom = sem_geom(Basis,Basisd,xc,N,Nd,nel,dxm1,dxtm1);

#L  = AssembleMatrix(c0,Geom.cnv,ν,Geom.wlp,lx1,nel);

#L  = AssembleMatrixLesshafft(c0,Geom.cnv,Geom.wlp,Geom.xm1,Geom.bm1,Basis,lx1,nel);

npts  = lx1*nel

ndof, glnum = Sem_Global_Num(Geom.xm1)

Q,QT   = Sem_QQT(glnum)
vmult  = sum(Q*QT,dims=2)
vimult = 1.0 ./vmult

xall  = vimult.*(Q*QT*Geom.xm1[:]);

c0 = 0.0e-10;

L,B,OP,Conv,Src,Lap,Fd = AssembleMatrixLesshafft2(c0,Geom.cnv,Geom.wlp,Geom.xm1,Geom.bm1,Basis,lx1,nel);

ifglobal = true

if ifglobal
  Cg    = QT*Conv*Q    # Global Convection matrix
  Lg    = QT*Lap*Q     # Global Laplacian matrix
  Sg    = QT*Src*Q     # Global Src matrix
  Fg    = QT*Fd*Q      # Global Feedback matrix
  Bg    = QT*B         # Global Mass vector
  Big   = 1.0./Bg      # Global inverse Mass vector
  
  OPg   = QT*(L)*Q./Bg
end

#OPg  = similar(Cg)

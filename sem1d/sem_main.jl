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


#L  = AssembleMatrix(c0,Geom.cnv,ν,Geom.wlp,lx1,nel);

#L  = AssembleMatrixLesshafft(c0,Geom.cnv,Geom.wlp,Geom.xm1,Geom.bm1,Basis,lx1,nel);

npts  = lx1*nel

ndof, glnum = Sem_Global_Num(Geom.xm1)

Q,QT   = Sem_QQT(glnum)
vmult  = sum(Q*QT,dims=2)
vimult = 1.0 ./vmult

xall  = vimult.*(Q*QT*Geom.xm1[:]);

c0 = 0.0;

L,B,OP,Conv,Src,Lap,Fd = AssembleMatrixLesshafft2(c0,Geom.cnv,Geom.wlp,Geom.xm1,Geom.bm1,Basis,lx1,nel);


println("Main interface for 1D SEM")

using PolynomialBases
using LinearAlgebra
using SparseArrays
using Printf

# Include the function files

include("$SRC/sem_geom.jl")
include("$SRC/AssembleMatrix.jl")
include("$SRC/Sem_QQT.jl")

# Load Parameters
#include("custom_params.jl")

# Generate the geomerty dependent matrices
Scale = 1.0
linf  = false
rinf  = true

#Geom  = sem_geom(Basis,Basisd,xc,N,Nd,nel,dxm1,dxtm1,linf,rinf,Scale,prec);
#Geom  = sem_geom(Basis,Basisd,xc,N,Nd,nel,dxm1,dxtm1,prec);

#Basis = SpectralBases.LaguerreSpectral(N)
Geom  = sem_geom_laguerre(Basis,M2N,N2M,prec)

npts  = lx1*nel
ndof  = npts

#ndof, glnum = Sem_Global_Num(Geom.xm1,prec,ifperiodic)

Q      = diagm(ones(VT,lx1))
QT     = Q'
vmult  = sum(Q*QT,dims=2)
vmult  = vmult[:]
vimult = one./vmult

vmultg  = sum(QT*Q,dims=2)
vmultg  = vmultg[:]
vimultg = one./vmultg


@printf("Done.")



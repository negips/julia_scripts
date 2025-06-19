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
Geom = sem_geom(Basis,Basisd,xc,N,Nd,nel,dxm1,dxtm1,prec);

npts  = lx1*nel
#ifperiodic = true
ndof, glnum = Sem_Global_Num(Geom.xm1,prec,ifperiodic)

Q,QT   = Sem_QQT(glnum,prec)
vmult  = sum(Q*QT,dims=2)
vmult  = vmult[:]
vimult = one./vmult

vmultg  = sum(QT*Q,dims=2)
vmultg  = vmultg[:]
vimultg = one./vmultg

#ifsparse  = true
L,B,OP,Conv,Src,Lap = AssembleMatrixCRD(U,γ,Geom.cnv,Geom.wlp,Geom.xm1,Geom.bm1,Basis,lx1,nel,prec,ifsparse);

# Build Dealiased Mass Matrix
if (ifsparse)
  Md    = spzeros(VT,npts,npts)
else
  Md    = zeros(VT,npts,npts)
end

for i in 1:nel
  j1 = (i-1)*lx1 + 1;
  j2 = i*lx1;
  # Dealiased Mass matrix
  Md[j1:j2,j1:j2] = (Geom.intpm1d')*diagm(Geom.bm1d[:,i])*Geom.intpm1d
end

# Build Filter Matrix
Fil,OPf = BuildFilter(M2N,N2M,lx1,nel,prec,ifsparse)

# Interpolation Operator
Intp    = AssembleInterpolation(Geom.intpm1d,nel,ifsparse)

#ifglobal = true

ndof2 = ndof + ndof
# Build Matrix for 2 fields together
if (ifsparse)
  L2g    = spzeros(VT,ndof2,ndof2)
  Z2g    = spzeros(VT,ndof,ndof)
else
  L2g    = zeros(VT,ndof2,ndof2)
  Z2g    = zeros(VT,ndof,ndof)
end

if ifglobal
  Cg    = QT*Conv*Q    # Global Convection matrix
  Lg    = QT*Lap*Q     # Global Laplacian matrix
  Sg    = QT*Src*Q     # Global Src matrix
  Filg  = QT*Fil*Q     # Global Filter matrix
  Intpg = Intp*Q       # Global Interpolation matrix

  Bg    = QT*B
  Big   = one./Bg      # Global inverse Mass vector
  Mdg   = QT*Md*Q      # Global Dialiased Weight Matrix for inner products 
 
  OPg   = QT*(L)*Q./Bg

  # Operator for 2 fields together.
  L2g   = [Z2g Z2g;γall[1]*Lg Z2g]
  B2g   = [Bg; Bg]
  @printf("Direct Global Matrices Built\n")
  
end



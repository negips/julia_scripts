println("Main interface for 1D SEM")

using PolynomialBases
using LinearAlgebra
using SparseArrays
# using UnicodePlots
#using Plots

# Include the function files
include("sem_geom.jl")

include("AssembleMatrix.jl")

include("AssembleMatrixLesshafft.jl")

include("sem_init_ref.jl")

include("Sem_QQT.jl")

Geom = sem_geom(Basis,Basisd,xc,N,Nd,nel,dxm1,dxtm1,prec);

npts  = lx1*nel
ndof, glnum = Sem_Global_Num(Geom.xm1,prec)

Q,QT   = Sem_QQT(glnum,prec)
vmult  = sum(Q*QT,dims=2)
vmult  = vmult[:]
vimult = one./vmult

xall  = vimult.*(Q*QT*Geom.xm1[:]);

if prec == BigFloat
  U   = BigFloat(6.0)
  γ   = BigFloat(1.0) - BigFloat(1.0)*im
  c0  = BigFloat(1.0e+02)
else
  U   = 6.0
  γ   = 1.0 - 1.0*im
  c0  = 1.0e+02
end  

# So much faster with Sparse Arrays
# Numerical error also seems much better behaved
ifsparse = true
if (ifsparse)
  L,B,OP,Conv,Src,Lap,Fd = AssembleMatrixLesshafftSparse(U,γ,c0,Geom.cnv,Geom.wlp,Geom.xm1,Geom.bm1,Basis,lx1,nel,prec);
else
  L,B,OP,Conv,Src,Lap,Fd = AssembleMatrixLesshafft2(U,γ,c0,Geom.cnv,Geom.wlp,Geom.xm1,Geom.bm1,Basis,lx1,nel,prec);
end  

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


#ifglobal = false

if ifglobal
  Cg    = QT*Conv*Q    # Global Convection matrix
  Lg    = QT*Lap*Q     # Global Laplacian matrix
  Sg    = QT*Src*Q     # Global Src matrix
  Fg    = QT*Fd*Q      # Global Feedback matrix
  Bg    = QT*B         # Global Mass vector
  Big   = one./Bg      # Global inverse Mass vector
  Mdg   = QT*Md*Q      # Global Dialiased Weight Matrix for inner products 
  
  OPg   = QT*(L)*Q./Bg
  
  println("Global Matrices Built")
end









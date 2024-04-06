println("Main interface for 1D SEM")


include("sem_set_geom.jl")

#ifsparse  = true
L,B,OP,Conv,Src,Lap = AssembleMatrixCRD(U,γ,Geom.cnv,Geom.wlp,Geom.xm1,Geom.bm1,Basis,lx1,nel,prec,ifsparse);

# Build Dealiased Mass Matrix
if (ifsparse)
  Md    = spzeros(VT,npts,npts)
else
  Md    = zeros(VT,npts,npts)
end

# Build Filter Matrix
Fil,OPf = BuildFilter(M2N,N2M,lx1,nel,prec,ifsparse)

#ifglobal = true

if ifglobal
  Cg    = QT*Conv*Q    # Global Convection matrix
  Lg    = QT*Lap*Q     # Global Laplacian matrix
  Sg    = QT*Src*Q     # Global Src matrix
  Filg  = QT*Fil*Q     # Global Filter matrix

#  Fg    = QT*Fd*Q      # Global Feedback matrix
#  BgM   = QT*diagm(B)*Q         # Global Mass vector
  Bg    = QT*B
  Big   = one./Bg      # Global inverse Mass vector
  Mdg   = QT*Md*Q      # Global Dialiased Weight Matrix for inner products 
  
  OPg   = QT*(L)*Q./Bg
  @printf("Direct Global Matrices Built\n")

end



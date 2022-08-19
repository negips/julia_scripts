println("Testing Fast Diagonalization Method")

using PolynomialBases
using LinearAlgebra
using PyPlot

include("sem_geom.jl")
#using PyPlot,PyCall

# Include the function files

#VT = Float64
#VT = ComplexF64
#setprecision(128)
#prec = BigFloat
prec = Float64

VT  = Complex{prec}

one = prec(1.0)
zro = prec(0.0)

# define nodal bases
N           = 7 ;                               # polynomial degree
lx1         = N+1;                              # No of points
Basis       = LobattoLegendre(N, prec)          # Polynomial Basis

Nd          = Int64(floor(N*3)+1)             # polynomial degree
lx1d        = Nd+1;                             # No of points
Basisd      = LobattoLegendre(Nd, prec)         # Polynomial Basis

N2          = N-2;                              # polynomial degree (Mesh 2)
lx2         = N2+1;                             # No of points
Basis2      = GaussLegendre(N2, prec)           # Polynomial Basis

xs          = prec(-1.)                          # Domain start
xe          = prec(1.0)                         # Domain end

ys          = prec(-1.)                          # Domain start
ye          = prec(1.0)                         # Domain end

nel         = 1                                 # No of elements
nnodes      = nel+1;                            # No of nodes
xc          = range(xs,stop=xe,length=nnodes)   # Element coordinates
yc          = range(ys,stop=ye,length=nnodes)   # Element coordinates

zgm1        = Basis.nodes;                      # Reference GLL Points
wzm1        = Basis.weights;                    # Reference integration weights

zgm1d       = Basisd.nodes;                     # Reference GLL Points
wzm1d       = Basisd.weights;                   # Reference integration weights

zgm2        = Basis2.nodes;                     # Reference GL Points
wzm2        = Basis2.weights;                   # Reference integration weights (GL)

# Reference Matrices
dxm1        = Basis.D;                          # Derivative Matrix
dxtm1       = Basis.D';                         # Derivative transpose

dxm2        = Basis2.D;                         # Derivative Matrix
dxtm2       = Basis2.D';                        # Derivative transpose

println("Reference Matrices Initialized")

Geomx  = sem_geom(Basis,Basisd,xc,N,Nd,nel,dxm1,dxtm1,prec)
Geom2x = sem_geom(Basis2,Basis,xc,N2,N,nel,dxm2,dxtm2,prec)
Geom12x = sem_geom(Basis,Basis2,xc,N,N2,nel,dxm1,dxtm1,prec)
Geomy  = sem_geom(Basis,Basisd,yc,N,Nd,nel,dxm1,dxtm1,prec)
Geom2y = sem_geom(Basis2,Basis,yc,N2,N,nel,dxm2,dxtm2,prec)


XM1 = zeros(prec,lx1,lx1,nel)
YM1 = zeros(prec,lx1,lx1,nel)
LAPx = zeros(prec,lx1,lx1,nel)
BMx  = zeros(prec,lx1,lx1,nel)

LAPy = zeros(prec,lx1,lx1,nel)
BMy  = zeros(prec,lx1,lx1,nel)

As   = zeros(prec,lx1*lx1,lx1*lx1,nel)
Bs   = zeros(prec,lx1*lx1,lx1*lx1,nel)
λxy  = zeros(prec,lx1*lx1,lx1*lx1,nel)
λinv = zeros(prec,lx1*lx1,lx1*lx1,nel)
Sx   = zeros(prec,lx1,lx1,nel)
STx  = zeros(prec,lx1,lx1,nel)
λx   = zeros(prec,lx1,nel)
Sy   = zeros(prec,lx1,lx1,nel)
STy  = zeros(prec,lx1,lx1,nel)
λy   = zeros(prec,lx1,nel)


one = ones(prec,1,lx1)
oneMx = Matrix{prec}(I,lx1,lx1)
oneMy = Matrix{prec}(I,lx1,lx1)

x0    = zeros(prec,lx1,lx1,nel)

for e=1:nel
  x1 = Geomx.xm1[:,e]
  XM1[:,:,e] = kron(x1,one)

  y1 = Geomy.xm1[:,e]'
  YM1[:,:,e] = kron(y1,one)

  bm = Geomx.bm1[:,1];
  BM = Diagonal(bm);
  LAPx[:,:,e] = Geomx.gradx[:,:,e]'*BM*Geomx.gradx[:,:,e]
  BMx[:,:,e]  = BM

  bm = Geomy.bm1[:,1];
  BM = Diagonal(bm);
  LAPy[:,:,e] = Geomy.gradx[:,:,e]'*BM*Geomy.gradx[:,:,e]
  BMy[:,:,e]  = BM

  
  As[:,:,e]   = kron(BMy[:,:,e],LAPx[:,:,e]) + kron(LAPy[:,:,e],BMx[:,:,e])
  Bs[:,:,e]   = kron(BMy[:,:,e],BMx[:,:,e])

  λx[:,e],Sx[:,:,e] = eigen(LAPx[:,:,e],BMx[:,:,e])
  STx[:,:,e] = Sx[:,:,e]'

  λy[:,e],Sy[:,:,e] = eigen(LAPy[:,:,e],BMy[:,:,e])
  STy[:,:,e] = Sy[:,:,e]'

  λxy[:,:,e] = kron(oneMy,Diagonal(λx[:,e])) + kron(Diagonal(λy[:,e]),oneMx)

  x0[:,:,e]  = kron(Geomy.bm1[:,e]',Geomx.bm1[:,e])
end  

for i in 1:(lx1*lx1)
  s = λxy[i,i]
  if abs(s)<1.0e-12
    s = 0.0
  else
    s = 1.0/s
  end

  λinv[i,i] = s
end

λbig,Sbig = eigen(As[:,:,1],Bs[:,:,1])


# This is what Nek does:
dgl = zeros(prec,lx2,lx1)
for i in 1:lx2
  z = zgm2[i]
  for j in 1:lx1
    v = zeros(prec,lx1)
    v[j] = 1.0
    dgl[i,j] = derivative_at(z,v,zgm1,Basis.baryweights)
  end
end  

bm2 = Diagonal(wzm2)
bdgl = bm2*dgl

lap2 = dgl'*bdgl 

lapnek = zeros(prec,lx1,lx1)
b0 = 0
b1 = (lx1-1) - 0

n = lx1-1
for j=1:n-1
  for i=1:n-1
    s = 0.0
    for k=b0:b1
      bb = (1.0/wzm1[k+1])
      s1 = bdgl[i,k+1]*bb*bdgl[j,k+1]
#      println("$i $j $k $bb $s1")
      s = s + s1 
    end  
    lapnek[i+1,j+1] = s
  end
end
lapnek[1,1]       = 1.0
lapnek[lx1,lx1]   = 1.0

#Interpolation operator to Mesh2
intpm12 = zeros(VT,lx1,lx1);

jgl  = interpolation_matrix(zgm2,Basis.nodes,Basis.baryweights);
bjgl = bm2*jgl  

massnek = zeros(prec,lx1,lx1)
b0 = 0
b1 = (lx1-1) - 0

n = lx1-1
for j=1:n-1
  for i=1:n-1
    s = 0.0
    for k=b0:b1
      bb = (1.0/wzm1[k+1])
      s1 = bjgl[i,k+1]*bb*bjgl[j,k+1]
#      println("$i $j $k $bb $s1")
      s = s + s1 
    end  
    massnek[i+1,j+1] = s
  end
end
massnek[1,1]       = 1.0
massnek[lx1,lx1]   = 1.0


# Nek Declaration
#g(0:lx1-1,0:lx1-1)
#jgl(1:lx2,0:lx1-1)
#bm(0:lx1-1)
#
#call rzero(g,(n+1)*(n+1))
#do j=1,n-1
#   do i=1,n-1
#      do k=b0,b1
#         g(i,j) = g(i,j) + gmm*jgl(i,k)*bm(k)*jgl(j,k)
#      enddo
#   enddo
#enddo


# Solving the linear System
α     = 2.0
x     = XM1[:,1];
x12   = copy(x);
x12[2:lx1-1] = Geom2x.xm1;
f     = -Geomx.bm1[:,1].*((α*π)^2).*sin.(α*π*x);
A     = LAPx[:,:,1];
B     = BMx[:,:,1];
# Set Boundary conditions:
d      = 1.0
A[1,:] = zeros(prec,1,lx1)
A[1,1] = 1.0
B[1,:] = zeros(prec,1,lx1)
B[1,1] = 1.0
f[1]   = d

A[lx1,:]    = zeros(prec,1,lx1)
A[lx1,lx1]  = 1.0
B[lx1,:]    = zeros(prec,1,lx1)
B[lx1,lx1]  = 1.0
f[lx1]      = d

u           = inv(A)*f
unek        = inv(lapnek)*f

plot(x,u)
plot(x12, unek)
plot(x,d .- sin.(α*π*x),linestyle="--",linewidth=2)

println("Done")















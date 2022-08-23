println("Testing Consistent Poisson Operator")

using PolynomialBases
using FastGaussQuadrature
using LinearAlgebra
using PyPlot
using IterativeSolvers

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

close("all")

# define nodal bases
N           = 7 ;                               # polynomial degree
lx1         = N+1;                              # No of points
Basis       = LobattoLegendre(N, prec)          # Polynomial Basis

Nd          = 6                                 # polynomial degree
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

dxm12 = zeros(prec,lx2,lx1)
for i in 1:lx2
  zi = zgm2[i]
  for j in 1:lx1
    v = zeros(prec,lx1)
    v[j] = 1.0
    dxm12[i,j] = derivative_at(zi,v,zgm1,Basis.baryweights)
  end
end  
dxtm12 = transpose(dxm12)


println("Reference Matrices Initialized")

G1x  = sem_geom(Basis,Basis2,xc,N,N2,nel,dxm1,dxtm1,prec)
G2x  = sem_geom(Basis2,Basis,xc,N2,N,nel,dxm2,dxtm2,prec)
G1y  = sem_geom(Basis,Basis2,yc,N,N2,nel,dxm1,dxtm1,prec)
G2y  = sem_geom(Basis2,Basis,yc,N2,N,nel,dxm2,dxtm2,prec)

G1d  = sem_geom(Basis,Basisd,xc,N,Nd,nel,dxm1,dxtm1,prec)



XM1   = zeros(prec,lx1,lx1,nel)
LAPx  = zeros(prec,lx1,lx1,nel)
BMx   = zeros(prec,lx1,lx1,nel)

XM2   = zeros(prec,lx2,lx2,nel)
BM2x  = zeros(prec,lx2,lx2,nel)

CPoisson = zeros(prec,lx2,lx2,nel)


one = ones(prec,1,lx1)
oneMx = Matrix{prec}(I,lx1,lx1)
oneMy = Matrix{prec}(I,lx1,lx1)

oneM2 = ones(prec,1,lx2)



for e=1:nel

  x1 = G1x.xm1[:,e]
  XM1[:,:,e] = kron(x1,one)

  bm = G1x.bm1[:,e];
  BM = Diagonal(bm);
  LAPx[:,:,e] = G1x.gradx[:,:,e]'*BM*G1x.gradx[:,:,e]
  BMx[:,:,e]  = BM

  x2  = G2x.xm1[:,e]
  XM2[:,:,e] = kron(x2,oneM2)
  bm2 = G2x.bm1[:,e]
  BM2 = Diagonal(bm2)
  BM2x[:,:,e]  = BM2

  rxm12 = G1d.intpm1d[:,:,e]*G1x.rxm1[:,e]

  gradm12 = Diagonal(rxm12)*dxm12
  gradtm12 = transpose(gradm12)
  bgradtm12 = gradtm12*BM2

  w1inv   = Diagonal(1.0./bm)
  CPoisson[:,:,e] = BM2*gradm12*w1inv*bgradtm12

end  

Fcp = eigen(CPoisson[:,:,1],BM2x[:,:,1])

x  = G2x.xm1[:,1]

f  = (1.0e-2)*x
bf = G2x.bm1[:,1].*f

A = CPoisson[:,:,1]

sol = gmres(A,bf)

Nf = 50
xf = Vector(range(xs,xe,length=Nf))
jf = interpolation_matrix(xf,Basis2.nodes,Basis2.baryweights);
solf = jf*sol

#plot(xf,solf)
#plot(x,sol)


# 2D Operator

p = (1.0e-2)*XM2[:,:,1]




#plot(xf,soln,linestyle="--",linewidth=2,label="Analytical")
#legend(fontsize=8)



println("Done")















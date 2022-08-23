println("Testing Consistent Poisson Operator")

using PolynomialBases
using FastGaussQuadrature
using LinearAlgebra
using PyPlot
using IterativeSolvers

include("sem_geom.jl")
include("sem_geom2d.jl")
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


G2D_M1  = sem_geom2d(Basis,Basis,xc,yc,nel,prec)
G2D_M2  = sem_geom2d(Basis2,Basis2,xc,yc,nel,prec)


dxm12 = zeros(prec,lx2,lx1)
for j in 1:lx1
  v = zeros(prec,lx1)
  v[j] = 1.0
  dxm12[:,j] = derivative_at(zgm2,v,Basis.nodes,Basis.baryweights)
end
dxtm12 = transpose(dxm12)

intm12 = zeros(prec,lx2,lx1)
intm12 = interpolation_matrix(Basis2.nodes,Basis.nodes,Basis.baryweights)

e = 1
p = (1.0e-2).*G2D_M2.x[:,:,e]

cdtp  = dxtm12*(G2D_M2.drdx[:,:,e].*p)



Nf = 50
xf = Vector(range(xs,xe,length=Nf))
jf = interpolation_matrix(xf,Basis2.nodes,Basis2.baryweights);
#solf = jf*sol

#plot(xf,solf)
#plot(x,sol)






#plot(xf,soln,linestyle="--",linewidth=2,label="Analytical")
#legend(fontsize=8)



println("Done")















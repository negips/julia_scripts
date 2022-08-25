println("Testing Consistent Poisson Operator")

using PolynomialBases
using FastGaussQuadrature
using LinearAlgebra
using PyPlot
using IterativeSolvers

include("sem_geom.jl")
include("sem_geom2d.jl")
include("ConsistentPoisson2D.jl")
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
N           = 20 ;                              # polynomial degree
lx1         = N+1;                              # No of points
Basis       = LobattoLegendre(N, prec)          # Polynomial Basis

N2          = N-2;                              # polynomial degree (Mesh 2)
lx2         = N2+1;                             # No of points
Basis2      = GaussLegendre(N2, prec)           # Polynomial Basis

N3          = N+10 ;                              # polynomial degree
lx3         = N3+1;                              # No of points
Basis3      = GaussLegendre(N3, prec)            # Polynomial Basis

N4          = floor(Int64,1.5*N) ;                              # polynomial degree
lx4         = N4+1;                              # No of points
Basis4      = GaussLegendre(N4, prec)            # Polynomial Basis

N5          = 2*N ;                              # polynomial degree
lx5         = N5+1;                              # No of points
Basis5      = GaussLegendre(N5, prec)            # Polynomial Basis


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

zgm2        = Basis2.nodes;                     # Reference GL Points
wzm2        = Basis2.weights;                   # Reference integration weights (GL)

zgm3        = Basis3.nodes;                     # Reference GLL Points
wzm3        = Basis3.weights;                   # Reference integration weights


# Reference Matrices
dxm1        = Basis.D;                          # Derivative Matrix
dxtm1       = Basis.D';                         # Derivative transpose

dxm2        = Basis2.D;                         # Derivative Matrix
dxtm2       = Basis2.D';                        # Derivative transpose

dxm3        = Basis3.D;                         # Derivative Matrix
dxtm3       = Basis3.D';                        # Derivative transpose

println("Reference Matrices Initialized")

G1x  = sem_geom(Basis,Basis2,xc,N,N2,nel,dxm1,dxtm1,prec)
G2x  = sem_geom(Basis2,Basis,xc,N2,N,nel,dxm2,dxtm2,prec)
G1y  = sem_geom(Basis,Basis2,yc,N,N2,nel,dxm1,dxtm1,prec)
G2y  = sem_geom(Basis2,Basis,yc,N2,N,nel,dxm2,dxtm2,prec)


G2D_M1  = sem_geom2d(Basis,Basis,xc,yc,nel,prec)
G2D_M2  = sem_geom2d(Basis2,Basis2,xc,yc,nel,prec)
G2D_M3  = sem_geom2d(Basis3,Basis3,xc,yc,nel,prec)
G2D_M4  = sem_geom2d(Basis4,Basis4,xc,yc,nel,prec)
G2D_M5  = sem_geom2d(Basis5,Basis5,xc,yc,nel,prec)


#dxm12 = zeros(prec,lx2,lx1)
#for j in 1:lx1
#  v = zeros(prec,lx1)
#  v[j] = 1.0
#  dxm12[:,j] = derivative_at(Basis2.nodes,v,Basis.nodes,Basis.baryweights)
#end
#dxtm12 = transpose(dxm12)
#
#intm12 = zeros(prec,lx2,lx1)
#intm12 = interpolation_matrix(Basis2.nodes,Basis.nodes,Basis.baryweights)
#
#e     = 1
#p     = (1.0e-0).*G2D_M2.x[:,:,e].*G2D_M2.b[:,:,e]    # rhs is multiplied by BM2
#p_M1  = (1.0e-0).*G2D_M1.x[:,:,e]
#
#cdtp  = dxtm12*(G2D_M2.drdx[:,:,e].*G2D_M2.b[:,:,e].*p)*intm12
#
#cddtp = G2D_M2.b[:,:,e].*(dxm12*(cdtp)*(intm12'))
#
## as a kronecker product:
#wgradm12x  = kron(intm12',dxtm12)*Diagonal(G2D_M2.b[:])
#wgradm12y  = kron(dxtm12,intm12')*Diagonal(G2D_M2.b[:])
#g_bp       = wgradm12x*p[:]
#kroncdtp   = reshape(g_bp,lx1,lx1)        #  == cdtp (Works)
#
#binv      = Diagonal(1.0./G2D_M1.b[:])
#
#div12x   =  Diagonal(G2D_M2.b[:])*kron(intm12,dxm12)*binv
#div12y   =  Diagonal(G2D_M2.b[:])*kron(dxm12,intm12)*binv
#
#d_g_bp  = div12x*g_bp
#kroncddtp = reshape(d_g_bp,lx2,lx2)       # == cddtp (Works)
#
#
#kroncddt  = div12x*wgradm12x + div12y*wgradm12y
#kronmass  = Matrix(Diagonal(G2D_M2.b[:]))

gradnodes   = Basis2.nodes
gradweights = Basis2.weights
divnodes    = Basis2.nodes
divweights  = Basis2.weights

binv = 1.0./G2D_M1.b
CPoisson2D_nek,Mass_nek = ConsistentPoisson2D(Basis,Basis2,binv,gradnodes,gradweights,divnodes,divweights,prec)

#CPoisson2D_nek = kroncddt

Fcp_nek   = eigen(CPoisson2D_nek,Matrix(Mass_nek))

## Consistent Integration (2D)
#--------------------------------------------------

gradnodes   = Basis3.nodes
gradweights = Basis3.weights
divnodes    = Basis3.nodes
divweights  = Basis3.weights
binv = 1.0./G2D_M1.b
CPoisson2D_132,Mass_132 = ConsistentPoisson2D(Basis,Basis2,binv,gradnodes,gradweights,divnodes,divweights,prec)

Fcp_132     = eigen(CPoisson2D_132,Matrix(Mass_132))

# Consistent Poisson with Dealiased Integration (2D) only for product
#--------------------------------------------------

gradnodes   = Basis3.nodes
gradweights = Basis3.weights
divnodes    = Basis3.nodes
divweights  = Basis3.weights
dealnodes   = Basis5.nodes
dealweights = Basis5.weights

CPoisson2D_13d2,Mass_13d2 = ConsistentPoisson2D_Dealias(Basis,Basis2,binv,gradnodes,gradweights,divnodes,divweights,dealnodes,dealweights,prec)
Fcp_13d2  = eigen(CPoisson2D_13d2,Matrix(Mass_13d2))


#-------------------------------------------------- 


plot(Fcp_nek.values,linestyle="none",marker="o",label="nek")
plot(Fcp_132.values,linestyle="none",marker="o",label="132")
#plot(Fcp_1342.values,linestyle="none",marker="o",label="D. 1342")
#plot(Fcp_1442.values,linestyle="none",marker="o",label="D. 1442")
plot(Fcp_13d2.values,linestyle="none",marker="o",label="D. 13d2")

legend(fontsize=12)


f          = 1.0*G2D_M2.x[:]
f2         = (1.0e-3)*ones(prec,lx2*lx2)

bf         = G2D_M2.b[:].*f

sol        = gmres(CPoisson2D_132,bf)
solm       = reshape(sol,lx2,lx2)

solnek     = gmres(CPoisson2D_nek,bf)
solmnek    = reshape(solnek,lx2,lx2)

sold3      = gmres(CPoisson2D_13d2,bf)
sold3m     = reshape(sold3,lx2,lx2)

Nf = 50
zxf = Vector(range(xs,xe,length=Nf))
jf = interpolation_matrix(zxf,Basis2.nodes,Basis2.baryweights);

xf = jf*G2D_M2.x[:,:,1]*(jf');
yf = jf*G2D_M2.y[:,:,1]*(jf');

sfm = jf*solm*(jf');
sfnek = jf*solmnek*(jf');
sfd3 = jf*sold3m*(jf');

#surf(G2D_M2.x[:,:,1],G2D_M2.y[:,:,1],(soldm .- solmnek))
#surf(G2D_M2.x[:,:,1],G2D_M2.y[:,:,1],soldm)
#surf(xf,yf,(sfd3 .- sfnek))

println("Done")















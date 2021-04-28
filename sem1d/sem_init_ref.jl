println("Intializing SEM")

#using PolynomialBases
#using PyPlot,PyCall

# Include the function files
#include("sem_geom.jl")

#close("all")

AT = Float64

# define nodal bases
N           = 16 ;                               # polynomial degree
lx1         = N+1;                              # No of points
Basis       = LobattoLegendre(N)                # Polynomial Basis

Nd          = 20 ;                              # polynomial degree
lx1d        = Nd+1;                             # No of points
Basisd      = LobattoLegendre(Nd)               # Polynomial Basis

#basis2 = GaussLegendre(N)

xs          = 0.                                # Domain start
xe          = 40.                                # Domain end
nel         = 20                                 # No of elements
nnodes      = nel+1;                            # No of nodes
xc          = range(xs,stop=xe,length=nnodes);  # Element coordinates

zgm1        = Basis.nodes;                      # Reference GLL Points
wzm1        = Basis.weights;                    # Reference integration weights

zgm1d       = Basisd.nodes;                     # Reference GLL Points
wzm1d       = Basisd.weights;                   # Reference integration weights


# Reference Matrices
dxm1        = Basis.D;                          # Derivative Matrix
dxtm1       = Basis.D';                         # Derivative transpose


#xm1, xrm1, rxm1, jacm1, jacmi, bm1, gradx, intpm1d, gradxd, bm1d, wlp, dvdx = sem_geom(Basis,Basisd,xc,N,Nd,nel,dxm1)
#
Geom = sem_geom(Basis,Basisd,xc,N,Nd,nel,dxm1,dxtm1);

println("Initialization done")

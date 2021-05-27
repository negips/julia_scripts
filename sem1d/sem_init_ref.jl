println("Intializing SEM")

#using PolynomialBases
#using PyPlot,PyCall

# Include the function files
#include("sem_geom.jl")

#close("all")

AT = Float64

# define nodal bases
N           = 8 ;                               # polynomial degree
lx1         = N+1;                              # No of points
Basis       = LobattoLegendre(N)                # Polynomial Basis

Nd          = Int64(floor(N*1.5)+1) ;                               # polynomial degree
lx1d        = Nd+1;                             # No of points
Basisd      = LobattoLegendre(Nd)               # Polynomial Basis

#basis2 = GaussLegendre(N)

xs          = 0.                                # Domain start
xe          = 50                              # Domain end
nel         = 30                                # No of elements
nnodes      = nel+1;                            # No of nodes
xc          = range(xs,stop=xe,length=nnodes);  # Element coordinates

zgm1        = Basis.nodes;                      # Reference GLL Points
wzm1        = Basis.weights;                    # Reference integration weights

zgm1d       = Basisd.nodes;                     # Reference GLL Points
wzm1d       = Basisd.weights;                   # Reference integration weights


# Reference Matrices
dxm1        = Basis.D;                          # Derivative Matrix
dxtm1       = Basis.D';                         # Derivative transpose

# Clean up matrices
eps10 = 10.0*eps(Float64)
for j in 1:lx1, i in 1:lx1
  if abs(dxm1[i,j])<eps10
    dxm1[i,j] = 0.
  end  
  if abs(dxtm1[i,j])<eps10
    dxtm1[i,j] = 0.
  end  
end


Geom = sem_geom(Basis,Basisd,xc,N,Nd,nel,dxm1,dxtm1);


println("Initialization done")

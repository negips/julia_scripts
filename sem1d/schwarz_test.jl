println("Additive Schwarz in SEM")

ENV["MPLBACKEND"]="qt5agg"

using PolynomialBases
using LinearAlgebra
using PyPlot
#using PyPlot,PyCall

# Include the function files

#VT = Float64
#VT = ComplexF64
#setprecision(90)
#prec = BigFloat
#prec = Double64
prec  = Float64

#VT    = Complex{prec}
VT    = Float64

one = prec(1.0)
zro = prec(0.0)

# define nodal bases
N           = 8 ;                               # polynomial degree
lx1         = N+1;                              # No of points
Basis       = LobattoLegendre(N, prec)          # Polynomial Basis

N2          = N-2 ;                             # polynomial degree
lx2         = N2+1;                             # No of points
Basis       = GaussLegendre(N, prec)            # Polynomial Basis


xs          = prec(0.)                          # Domain start
xe          = prec(20.0)                         # Domain end

nel         = 10                                # No of elements
nnodes      = nel+1;                            # No of nodes
xc          = range(xs,stop=xe,length=nnodes)   # Element coordinates

zgm1        = Basis.nodes;                      # Reference GLL Points
wzm1        = Basis.weights;                    # Reference integration weights


# Reference Matrices
dxm1        = Basis.D;                          # Derivative Matrix
dxtm1       = Basis.D';                         # Derivative transpose

println("Reference Matrices Initialized")

x = 1 .+ Basis.nodes
y = sin.(1.0*pi.*x)

close("all")

plot(x,y)

A    = Matrix(Diagonal(Basis.weights))
d1t  = Basis.weights[1].*dxm1[1,:];
d1   = d1t'

dnt  = Basis.weights[lx1].*dxm1[lx1,:];
dn   = dnt'


B   = [A  d1t dnt;
       d1  0.0 0.0;
       dn  0.0 0.0]

rhs = [Basis.weights.*y; 0.0; 0.0]

sol = inv(B)*rhs
z   = sol[1:lx1]

plot(x,z)

dxm1*z






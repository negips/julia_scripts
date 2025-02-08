println("Intializing SEM")

using PolynomialBases
#using PyPlot,PyCall

# Include the function files

#VT = Float64
#VT = ComplexF64
#setprecision(90)
#prec = BigFloat
#prec = Double64
prec  = Float64

VT  = Complex{prec}

one = prec(1.0)
zro = prec(0.0)

# define nodal bases
N           = 8  ;                              # polynomial degree
lx1         = N+1;                              # No of points
Basis       = LobattoLegendre(N, prec)          # Polynomial Basis

Nd          = Int64(floor(N*1.5)+1)             # polynomial degree
lx1d        = Nd+1;                             # No of points
Basisd      = LobattoLegendre(Nd, prec)         # Polynomial Basis

#basis2 = GaussLegendre(N)

xs          = prec(0.)                          # Domain start
xe          = prec(40.0)                        # Domain end

nel         = 41                                # No of elements
nnodes      = nel+1;                            # No of nodes
xc          = range(xs,stop=xe,length=nnodes)   # Element coordinates

zgm1        = Basis.nodes;                      # Reference GLL Points
wzm1        = Basis.weights;                    # Reference integration weights

zgm1d       = Basisd.nodes;                     # Reference GLL Points
wzm1d       = Basisd.weights;                   # Reference integration weights


# Reference Matrices
dxm1        = Basis.D;                          # Derivative Matrix
dxtm1       = Basis.D';                         # Derivative transpose


println("Reference Matrices Initialized")

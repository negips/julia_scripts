println("Intializing SEM")

using PolynomialBases
#using PyPlot,PyCall

# Include the function files

#VT = Float64
#VT = ComplexF64
setprecision(256)
prec = BigFloat
#prec = Float64

VT  = Complex{prec}

if prec == BigFloat
  one = BigFloat(1.0)
  zro = BigFloat(0.0)
else
  one = 1.0
  zro = 0.0
end  

# define nodal bases
N           = 11 ;                              # polynomial degree
lx1         = N+1;                              # No of points
Basis       = LobattoLegendre(N, prec)          # Polynomial Basis

Nd          = Int64(floor(N*1.5)+1)             # polynomial degree
lx1d        = Nd+1;                             # No of points
Basisd      = LobattoLegendre(Nd, prec)               # Polynomial Basis

#basis2 = GaussLegendre(N)

if (prec==BigFloat)
  xs          = BigFloat(0.)                    # Domain start
  xe          = BigFloat(40.0)                  # Domain end
else
  xs          = 0.                              # Domain start
  xe          = 40.0                            # Domain end
end  

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

# # Clean up matrices
# eps10 = 10.0*eps(Float64)
# for j in 1:lx1, i in 1:lx1
#   if abs(dxm1[i,j])<eps10
#     dxm1[i,j] = 0.
#   end  
#   if abs(dxtm1[i,j])<eps10
#     dxtm1[i,j] = 0.
#   end  
# end

println("Reference Matrices Initialized")

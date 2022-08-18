println("Testing The Schwarz Local Solves")

using PolynomialBases
using LinearAlgebra
using PyPlot
using IterativeSolvers

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
N           = 4 ;                               # polynomial degree
lx1         = N+1;                              # No of points
Basis       = LobattoLegendre(N, prec)          # Polynomial Basis

Nd          = Int64(floor(N*1.5)+1)             # polynomial degree
lx1d        = Nd+1;                             # No of points
Basisd      = LobattoLegendre(Nd, prec)         # Polynomial Basis

#basis2 = GaussLegendre(N)

xs          = prec(0.)                          # Domain start
xe          = prec(4.0)                         # Domain end

nel         = 2                                # No of elements
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

bm1 = Diagonal(wzm1)
wkLap = dxm1'*bm1*dxm1

bdry        = zeros(prec,lx1,lx1)
bdry[1,:]   = bm1[1,1]*dxm1[1,:]
bdry[lx1,:] = bm1[lx1,lx1]*dxm1[lx1,:]

wkLap[1,:]        = zeros(prec,1,lx1)
wkLap[1,1]        = 1.0

wkLap[lx1,:]      = zeros(prec,1,lx1)
wkLap[lx1,lx1]    = 1.0

#wkLap = zeros(prec,lx1,lx1)
#for i in 1:lx1
#  for j in 1:lx1
#    s = 0.0
#    for k in 1:lx1
#     s = s + dxm1[k,i]*wzm1[k]*dxm1[k,j]
#    end
#    wkLap[i,j] = s
#  end
#end
#
#Lap = zeros(prec,lx1,lx1)
#for i in 1:lx1
#  for j in 1:lx1
#    s = 0.0
#    for k in 1:lx1
#      s = s + wzm1[k]*dxm1[i,k]*dxm1[k,j]
#    end
#    Lap[i,j] = s
#  end
#end
#

x1 = zgm1 .+ 1.0
x2 = zgm1 .+ 3.0

v1 = ones(prec,lx1,1)
v2 = ones(prec,lx1,1)

#v1 = cos.(x1)
#v2 = cos.(x2)

#s1 = gmres(wkLap,v1,reltol=1.0e-6)

#s1 = pinv(wkLap)

#plot()








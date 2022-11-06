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
N           = 9 ;                              # polynomial degree
lx1         = N+1;                              # No of points
Basis       = LobattoLegendre(N, prec)          # Polynomial Basis

N2          = N-2;                              # polynomial degree (Mesh 2)
lx2         = N2+1;                             # No of points
Basis2      = GaussLegendre(N2, prec)           # Polynomial Basis

N3          = 3*N ;                              # polynomial degree
lx3         = N3+1;                              # No of points
Basis3      = LobattoLegendre(N3, prec)            # Polynomial Basis

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


## Nek Pseudo-Laplacian (2D)
#--------------------------------------------------

gradnodes   = Basis2.nodes
gradweights = Basis2.weights
divnodes    = Basis2.nodes
divweights  = Basis2.weights

binv = 1.0./G2D_M1.b
CPoisson2D_nek,Mass_nek = ConsistentPoisson2D(Basis,Basis2,binv,gradnodes,gradweights,divnodes,divweights,prec)


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
binv        = 1.0./G2D_M1.b

CPoisson2D_13d2,Mass_13d2 = ConsistentPoisson2D_Dealias(Basis,Basis2,binv,gradnodes,gradweights,divnodes,divweights,dealnodes,dealweights,prec)
Fcp_13d2  = eigen(CPoisson2D_13d2,Matrix(Mass_13d2))

## Nek Pseudo-Laplacian (2D) with variable density
#--------------------------------------------------

gradnodes   = Basis2.nodes
gradweights = Basis2.weights
divnodes    = Basis2.nodes
divweights  = Basis2.weights

x0          = -ones(prec,size(G2D_M1.x))
μ           = (1.0/0.7)
ρ           = exp.(-(μ.*(G2D_M1.x .- x0)).^2)
binv        = 1.0./(G2D_M1.b.*ρ)

CPoisson2D_nek_ρ,Mass_nek_ρ = ConsistentPoisson2D(Basis,Basis2,binv,gradnodes,gradweights,divnodes,divweights,prec)

Fcp_nek_ρ   = eigen(CPoisson2D_nek_ρ,Matrix(Mass_nek_ρ))

## Nek Pseudo-Laplacian (2D) with variable density
#--------------------------------------------------

gradnodes   = Basis2.nodes
gradweights = Basis2.weights
divnodes    = Basis2.nodes
divweights  = Basis2.weights

x0          = -ones(prec,size(G2D_M1.x))
μ           = (1.0/0.7)
ρ           = exp.(-(μ.*(G2D_M1.x .- x0)).^2)
binv        = 1.0./(G2D_M1.b)

CPoisson2D_nek_ρρ,Mass_nek_ρρ = ConsistentPoisson2D_ρ(Basis,Basis2,binv,ρ,gradnodes,gradweights,divnodes,divweights,prec)

Fcp_nek_ρρ   = eigen(CPoisson2D_nek_ρρ,Matrix(Mass_nek_ρρ))

## Consistent Integration (2D) with variable density
#--------------------------------------------------

gradnodes   = Basis3.nodes
gradweights = Basis3.weights
divnodes    = Basis3.nodes
divweights  = Basis3.weights

x0          = -ones(prec,size(G2D_M1.x))
μ           = (1.0/0.7)
ρ           = exp.(-(μ.*(G2D_M1.x .- x0)).^2)
binv        = 1.0./(G2D_M1.b.*ρ)
CPoisson2D_132_ρ,Mass_132_ρ = ConsistentPoisson2D(Basis,Basis2,binv,gradnodes,gradweights,divnodes,divweights,prec)

Fcp_132_ρ   = eigen(CPoisson2D_132_ρ,Matrix(Mass_132_ρ))

## Consistent Integration (2D) with variable density
#--------------------------------------------------

gradnodes   = Basis3.nodes
gradweights = Basis3.weights
divnodes    = Basis3.nodes
divweights  = Basis3.weights

x0          = -ones(prec,size(G2D_M1.x))
μ           = (1.0/0.7)
ρ           = exp.(-(μ.*(G2D_M1.x .- x0)).^2)
binv        = 1.0./(G2D_M1.b)

CPoisson2D_132_ρρ,Mass_132_ρρ = ConsistentPoisson2D_ρ(Basis,Basis2,binv,ρ,gradnodes,gradweights,divnodes,divweights,prec)

Fcp_132_ρρ   = eigen(CPoisson2D_132_ρρ,Matrix(Mass_132_ρρ))


# Consistent Poisson with Dealiased Integration (2D) only for product
#--------------------------------------------------

gradnodes   = Basis3.nodes
gradweights = Basis3.weights
divnodes    = Basis3.nodes
divweights  = Basis3.weights
dealnodes   = Basis5.nodes
dealweights = Basis5.weights
x0          = -ones(prec,size(G2D_M1.x))
μ           = (1.0/0.7)
ρ           = exp.(-(μ.*(G2D_M1.x .- x0)).^2)
binv        = 1.0./G2D_M1.b

CPoisson2D_13d2_ρ,Mass_13d2_ρ = ConsistentPoisson2D_Dealias_ρ(Basis,Basis2,binv,ρ,gradnodes,gradweights,divnodes,divweights,dealnodes,dealweights,prec)
Fcp_13d2_ρ    = eigen(CPoisson2D_13d2_ρ,Matrix(Mass_13d2_ρ))

# Projecting out the null-space
#-------------------------------------------------- 

p0  = ones(prec,length(G2D_M2.x))
vol = p0'*Mass_nek*p0
p0n = p0./sqrt(vol)
Pr0 = p0n*p0n'*Mass_nek       # Projector on to the null space
U   = I - Pr0                 # Projector orthogonal to the null space

CPoisson2D_nek0   = U*CPoisson2D_nek*U
Fcp_nek0   = eigen(CPoisson2D_nek0,Matrix(Mass_nek))


# Projecting out the null-space (without Mass Matrix)
#-------------------------------------------------- 

p01  = ones(prec,length(G2D_M2.x))
vol1 = p01'*p01 + 0.0
p01n = p01./sqrt(vol1)
Pr01 = p01n*p01n'             # Projector on to the null space
U1   = I - Pr01               # Projector orthogonal to the null space

CPoisson2D_nek1   = U1'*CPoisson2D_nek*U1
Fcp_nek1   = eigen(CPoisson2D_nek1,Matrix(Mass_nek))


# Projecting out the null-space of 132
#-------------------------------------------------- 

#p0  = ones(prec,length(G2D_M2.x))
#vol = p0'*Mass_nek*p0
#p0n = p0./sqrt(vol)
#Pr0 = p0n*p0n'*Mass_nek       # Projector on to the null space
#U   = I - Pr0                 # Projector orthogonal to the null space

CPoisson2D_1320   = U'*CPoisson2D_132*U
Fcp_1320   = eigen(CPoisson2D_1320,Matrix(Mass_132))


#-------------------------------------------------- 

<<<<<<< HEAD
#plot(Fcp_nek.values,linestyle="none",marker="o",label="nek")
#plot(Fcp_nek0.values,linestyle="none",marker="o",label="nek-0")
#plot(Fcp_nek1.values,linestyle="none",marker="o",label="nek-1")
plot(Fcp_nek_ρρ.values,linestyle="none",marker="o",label="nek-ρ")
#plot(Fcp_132_ρ.values,linestyle="none",marker="o",label="132-ρ")
plot(Fcp_132_ρρ.values,linestyle="none",marker="o",label="132-ρρ")
#plot(Fcp_13d2_ρ.values,linestyle="none",marker="o",label="13d2-ρ")
=======
plot(Fcp_nek.values,linestyle="none",marker="o",label="nek")
plot(Fcp_nek0.values,linestyle="none",marker="o",label="nek-B0")
plot(Fcp_nek1.values,linestyle="none",marker="o",label="nek-I0")
#plot(Fcp_1320.values,linestyle="none",marker="o",label="132-0")
>>>>>>> 429b6a0 (updating from workstation)
#plot(Fcp_132.values,linestyle="none",marker="o",label="132")
#plot(Fcp_1342.values,linestyle="none",marker="o",label="D. 1342")
#plot(Fcp_1442.values,linestyle="none",marker="o",label="D. 1442")
#plot(Fcp_13d2.values,linestyle="none",marker="o",label="D. 13d2")

legend(fontsize=12)


f          = 1.0*G2D_M2.x[:]
f2         = (1.0e-3)*ones(prec,lx2*lx2)

bf         = G2D_M2.b[:].*f

sol        = gmres(CPoisson2D_132,bf)
solm       = reshape(sol,lx2,lx2)

solnek     = gmres(CPoisson2D_nek,bf)
solmnek    = reshape(solnek,lx2,lx2)

solnek0    = gmres(CPoisson2D_nek0,bf)
solmnek0   = reshape(solnek0,lx2,lx2)

Nf = 50
zxf = Vector(range(xs,xe,length=Nf))
jf = interpolation_matrix(zxf,Basis2.nodes,Basis2.baryweights);

xf = jf*G2D_M2.x[:,:,1]*(jf');
yf = jf*G2D_M2.y[:,:,1]*(jf');

sfm    = jf*solm*(jf');
sfnek  = jf*solmnek*(jf');
sfnek0 = jf*solmnek0*(jf');

<<<<<<< HEAD
evnek = real.(reshape(Fcp_1320.vectors,lx2,lx2,lx2*lx2))

i = 1
evi   = evnek[:,:,i]; 

evf = jf*evi*(jf');
h2 = figure(num=2)
#surf(xf,yf,evf)


evnek_matrix = reshape(Fcp_nek_ρ.vectors,lx2,lx2,lx2*lx2)
evnew_matrix = reshape(Fcp_132_ρ.vectors,lx2,lx2,lx2*lx2)

evnek_matrix = reshape(Fcp_nek.vectors,lx2,lx2,lx2*lx2)
evnew_matrix = reshape(Fcp_132.vectors,lx2,lx2,lx2*lx2)
=======
#evnek = real.(reshape(Fcp_nek1.vectors,lx2,lx2,lx2*lx2))
#
#i = 1
#evi   = evnek[:,:,i]; 
#
#evf = jf*evi*(jf');
#h2 = figure(num=2)
#surf(xf,yf,evf)


evnek_matrix = real.(reshape(Fcp_nek0.vectors,lx2,lx2,lx2*lx2))
evnew_matrix = real.(reshape(Fcp_132.vectors,lx2,lx2,lx2*lx2))
>>>>>>> 429b6a0 (updating from workstation)

i = 1
evneki       = jf*evnek_matrix[:,:,i]*jf'
evnewi       = jf*evnew_matrix[:,:,i]*jf'

h2 = figure(num=2)
surf(xf,yf,evneki)

#surf(G2D_M2.x[:,:,1],G2D_M2.y[:,:,1],(soldm .- solmnek))
#surf(G2D_M2.x[:,:,1],G2D_M2.y[:,:,1],soldm)
#surf(xf,yf,sfnek)


println("Done")















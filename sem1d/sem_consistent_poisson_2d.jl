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
N           = 7 ;                              # polynomial degree
lx1         = N+1;                              # No of points
Basis       = LobattoLegendre(N, prec)          # Polynomial Basis

N2          = N-2;                              # polynomial degree (Mesh 2)
lx2         = N2+1;                             # No of points
Basis2      = GaussLegendre(N2, prec)           # Polynomial Basis

N3          = N+0 ;                              # polynomial degree
lx3         = N3+1;                              # No of points
Basis3      = LobattoLegendre(N3, prec)          # Polynomial Basis

N4          = N+8 ;                              # polynomial degree
lx4         = N4+1;                              # No of points
Basis4      = LobattoLegendre(N4, prec)          # Polynomial Basis

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


dxm12 = zeros(prec,lx2,lx1)
for j in 1:lx1
  v = zeros(prec,lx1)
  v[j] = 1.0
  dxm12[:,j] = derivative_at(Basis2.nodes,v,Basis.nodes,Basis.baryweights)
end
dxtm12 = transpose(dxm12)

intm12 = zeros(prec,lx2,lx1)
intm12 = interpolation_matrix(Basis2.nodes,Basis.nodes,Basis.baryweights)

e     = 1
p     = (1.0e-0).*G2D_M2.x[:,:,e].*G2D_M2.b[:,:,e]    # rhs is multiplied by BM2
p_M1  = (1.0e-0).*G2D_M1.x[:,:,e]

cdtp  = dxtm12*(G2D_M2.drdx[:,:,e].*G2D_M2.b[:,:,e].*p)*intm12

cddtp = G2D_M2.b[:,:,e].*(dxm12*(cdtp)*(intm12'))

# as a kronecker product:
wgradm12x  = kron(intm12',dxtm12)*Diagonal(G2D_M2.b[:])
wgradm12y  = kron(dxtm12,intm12')*Diagonal(G2D_M2.b[:])
g_bp       = wgradm12x*p[:]
kroncdtp   = reshape(g_bp,lx1,lx1)        #  == cdtp (Works)

binv      = Diagonal(1.0./G2D_M1.b[:])

div12x   =  Diagonal(G2D_M2.b[:])*kron(intm12,dxm12)*binv
div12y   =  Diagonal(G2D_M2.b[:])*kron(dxm12,intm12)*binv

d_g_bp  = div12x*g_bp
kroncddtp = reshape(d_g_bp,lx2,lx2)       # == cddtp (Works)


kroncddt  = div12x*wgradm12x + div12y*wgradm12y
kronmass  = Matrix(Diagonal(G2D_M2.b[:]))

CPoisson2D_nek = kroncddt

Fcp_nek   = eigen(kroncddt,kronmass)
#-------------------------------------------------- 

# Consistent Integration (1D)

dxm13 = zeros(prec,lx3,lx1)
for j in 1:lx1
  v = zeros(prec,lx1)
  v[j] = 1.0
  dxm13[:,j] = derivative_at(Basis3.nodes,v,Basis.nodes,Basis.baryweights)
end
dxtm13 = transpose(dxm13)

intm13 = zeros(prec,lx3,lx1)
intm13 = interpolation_matrix(Basis3.nodes,Basis.nodes,Basis.baryweights)

intm23 = zeros(prec,lx3,lx2)
intm23 = interpolation_matrix(Basis3.nodes,Basis2.nodes,Basis2.baryweights)
intm22 = Matrix{prec}(I,lx2,lx2)    # Just Identity.

wgradm231 = dxtm13*Diagonal(Basis3.weights)*intm23     # Basis2 -> Basis1 via Basis3 number of points
div132    = (intm23')*Diagonal(Basis3.weights)*dxm13   # Basis1 -> Basis2 via Basis3 number of points
CPoisson_new = div132*Diagonal((1.0./Basis.weights))*wgradm231

wgradm221 = dxtm12*Diagonal(Basis2.weights)*intm22     # Basis2 -> Basis1 via Basis2 number of points
div122    = (intm22')*Diagonal(Basis2.weights)*dxm12   # Basis1 -> Basis2 via Basis2 number of points
CPoisson_nek = div122*(Diagonal(1.0./Basis.weights))*wgradm221

# Consistent Integration (2D)
#--------------------------------------------------

# as a kronecker product:
intm23_2d  = kron(intm23,intm23);
wg2d_231_x = kron(intm13',dxtm13)*Diagonal(G2D_M3.b[:])*intm23_2d
wg2d_231_y = kron(dxtm13,intm13')*Diagonal(G2D_M3.b[:])*intm23_2d

wg2d_231 = [wg2d_231_x; wg2d_231_y]

binv      = 1.0./G2D_M1.b[:]

div2d_132_x  = (intm23_2d')*Diagonal(G2D_M3.b[:])*kron(intm13,dxm13)*Diagonal(binv)
div2d_132_y  = (intm23_2d')*Diagonal(G2D_M3.b[:])*kron(dxm13,intm13)*Diagonal(binv)

div2d_132 = [div2d_132_x div2d_132_y]


#CPoisson2D_new = div2d_132*wg2d_231
CPoisson2D_new = div2d_132_x*wg2d_231_x + div2d_132_y*wg2d_231_y

Fcp_new     = eigen(CPoisson2D_new,kronmass)

# Consistent Poisson with Dealiased Integration (2D)
#--------------------------------------------------

dxm14 = zeros(prec,lx4,lx1)
for j in 1:lx1
  v = zeros(prec,lx1)
  v[j] = 1.0
  dxm14[:,j] = derivative_at(Basis4.nodes,v,Basis.nodes,Basis.baryweights)
end
dxtm14 = transpose(dxm14)

intm14 = zeros(prec,lx4,lx1)
intm14 = interpolation_matrix(Basis4.nodes,Basis.nodes,Basis.baryweights)
intm14_2d  = kron(intm14,intm14);

intm24 = zeros(prec,lx2,lx4)
intm24 = interpolation_matrix(Basis4.nodes,Basis2.nodes,Basis2.baryweights)
intm24_2d  = kron(intm24,intm24);

intm44 = Matrix{prec}(I,lx4,lx4)    # Just Identity.

wg2d_231_x = kron(intm13',dxtm13)*Diagonal(G2D_M3.b[:])*intm23_2d
wg2d_231_y = kron(dxtm13,intm13')*Diagonal(G2D_M3.b[:])*intm23_2d

wg2d_231 = [wg2d_231_x; wg2d_231_y]

binv      = (1.0./G2D_M1.b[:])
binv4     = intm14_2d*binv

div2d_142_x  = (intm24_2d')*Diagonal(G2D_M4.b[:])*kron(intm44,Basis4.D)*Diagonal(binv4)*intm14_2d
div2d_142_y  = (intm24_2d')*Diagonal(G2D_M4.b[:])*kron(Basis4.D,intm44)*Diagonal(binv4)*intm14_2d

#div2d_142_x  = (intm24_2d')*Diagonal(G2D_M4.b[:])*kron(intm14,dxm14)*Diagonal(binv4)*intm14_2d
#div2d_142_y  = (intm24_2d')*Diagonal(G2D_M4.b[:])*kron(dxm14,intm14)*Diagonal(binv4)*intm14_2d

div2d_1342   = [div2d_142_x div2d_142_y]

CPoisson2D_1342 = div2d_142_x*wg2d_231_x + div2d_142_y*wg2d_231_y

Fcp_1342     = eigen(CPoisson2D_1342,kronmass)

#-------------------------------------------------- 

# Consistent Poisson with Dealiased Integration (2D) (both derivatives)
#--------------------------------------------------

wg2d_241_x = kron(intm14',dxtm14)*Diagonal(G2D_M4.b[:])*intm24_2d
wg2d_241_y = kron(dxtm14,intm14')*Diagonal(G2D_M4.b[:])*intm24_2d

wg2d_241 = [wg2d_241_x; wg2d_241_y]

binv      = (1.0./G2D_M1.b[:])
binv4     = intm14_2d*binv

div2d_142_x  = (intm24_2d')*Diagonal(G2D_M4.b[:])*kron(intm44,Basis4.D)*Diagonal(binv4)*intm14_2d
div2d_142_y  = (intm24_2d')*Diagonal(G2D_M4.b[:])*kron(Basis4.D,intm44)*Diagonal(binv4)*intm14_2d

#div2d_142_x  = (intm24_2d')*Diagonal(G2D_M4.b[:])*kron(intm14,dxm14)*Diagonal(binv4)*intm14_2d
#div2d_142_y  = (intm24_2d')*Diagonal(G2D_M4.b[:])*kron(dxm14,intm14)*Diagonal(binv4)*intm14_2d

div2d_1442   = [div2d_142_x div2d_142_y]

CPoisson2D_1442 = div2d_142_x*wg2d_241_x + div2d_142_y*wg2d_241_y

Fcp_1442     = eigen(CPoisson2D_1442,kronmass)

#-------------------------------------------------- 

# Consistent Poisson with Dealiased Integration (2D) only for product
#--------------------------------------------------

wg2d_231_x = kron(intm13',dxtm13)*Diagonal(G2D_M3.b[:])*intm23_2d
wg2d_231_y = kron(dxtm13,intm13')*Diagonal(G2D_M3.b[:])*intm23_2d

wg2d_231 = [wg2d_231_x; wg2d_231_y]

binv      = (1.0./G2D_M1.b[:])

deal      = Diagonal(1.0./G2D_M1.b[:])*(intm14_2d')*Diagonal(G2D_M4.b[:].*(intm14_2d*binv))*intm14_2d

div2d_13d2_x  = (intm23_2d')*Diagonal(G2D_M3.b[:])*kron(intm13,dxm13)*deal
div2d_13d2_y  = (intm23_2d')*Diagonal(G2D_M3.b[:])*kron(dxm13,intm13)*deal

#div2d_142_x  = (intm24_2d')*Diagonal(G2D_M4.b[:])*kron(intm14,dxm14)*Diagonal(binv4)*intm14_2d
#div2d_142_y  = (intm24_2d')*Diagonal(G2D_M4.b[:])*kron(dxm14,intm14)*Diagonal(binv4)*intm14_2d

div2d_13d2   = [div2d_13d2_x div2d_13d2_y]

CPoisson2D_13d2 = div2d_13d2_x*wg2d_231_x + div2d_13d2_y*wg2d_231_y

Fcp_13d2     = eigen(CPoisson2D_13d2,kronmass)

#-------------------------------------------------- 


plot(Fcp_nek.values,linestyle="none",marker="o",label="nek")
plot(Fcp_new.values,linestyle="none",marker="o",label="M3")
#plot(Fcp_1342.values,linestyle="none",marker="o",label="Dealias")
#plot(Fcp_1442.values,linestyle="none",marker="o",label="Dealias 2")
plot(Fcp_13d2.values,linestyle="none",marker="o",label="Dealias 3")

legend(fontsize=12)


f          = 1.0*G2D_M2.x[:]
f2         = (1.0e-3)*ones(prec,lx2*lx2)

bf         = G2D_M2.b[:].*f2

sol        = gmres(CPoisson2D_new,bf)
solm       = reshape(sol,lx2,lx2)

solnek     = gmres(kroncddt,bf)
solmnek    = reshape(solnek,lx2,lx2)

sold       = gmres(CPoisson2D_1342,bf)
soldm      = reshape(sold,lx2,lx2)


Nf = 50
zxf = Vector(range(xs,xe,length=Nf))
jf = interpolation_matrix(zxf,Basis2.nodes,Basis2.baryweights);

xf = jf*G2D_M2.x[:,:,1]*(jf');
yf = jf*G2D_M2.y[:,:,1]*(jf');

sfm = jf*solm*(jf');
sfnek = jf*solmnek*(jf');
sfd = jf*soldm*(jf');

#surf(G2D_M2.x[:,:,1],G2D_M2.y[:,:,1],(soldm .- solmnek))
#surf(G2D_M2.x[:,:,1],G2D_M2.y[:,:,1],soldm)
#surf(xf,yf,sfm)

println("Done")















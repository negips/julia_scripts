#!/bin/julia

println("Nearest Kronecker Product Approximation")

using LinearAlgebra
using PolynomialBases
using PyPlot
using Random
using Printf

#----------------------------------------------------------------------  
function ReArrangement(A,m1,m2,n1,n2)

  r,c = size(A)

  @assert m1*m2 == r && n1*n2 == c

  T  = eltype(A)
  RA = zeros(T,m1*n1,m2*n2)

  for j in 1:n1
    Rj  = zeros(T,m1,m2*n2)
    for i in 1:m1
      Aij    = A[(i-1)*m2+1:i*m2,(j-1)*n2+1:j*n2]
      Rjview = view(Rj,i,:)
      copyto!(Rjview,1,Aij,1,m2*n2)
    end
    RAview = view(RA,(j-1)*m1+1:j*m1,:)
    copyto!(RAview,1,Rj,1,m1*m2*n2)
  end

  return RA
end
#----------------------------------------------------------------------  



Nx    = 11
lx1   = Nx+1
Bx    = LobattoLegendre(Nx)
Vx    = ones(Float64,lx1)

Ny    = 11
ly1   = Ny+1
B     = LobattoLegendre(Ny)
Vy    = ones(Float64,ly1)

xs    = 0.0
xe    = 1.0
ys    = 1.0
ye    = 2.0

# X-direction
X     = zeros(Float64,lx1)
for i in 1:lx1
  X[i] = map_from_canonical(B.nodes[i],xs,xe,B)
end
X2D   = kron(Vy',X)

# Y-direction
Y     = zeros(Float64,ly1)
for i in 1:ly1
  Y[i] = map_from_canonical(B.nodes[i],ys,ye,B)
end
Y2D   = kron(Y',Vx)

rng   = MersenneTwister(1234)

Fld   = zeros(Float64,lx1,ly1)
for ci in CartesianIndices(Fld)
  x = X2D[ci]
  y = Y2D[ci]
  Fld[ci] = 0.0*rand(rng) + cos(π*x) + sin(π*y) + tan(x+y)
end

Fv = Fld[:]

FM = diagm(Fv)

m1 = ly1
n1 = ly1

m2 = lx1
n2 = lx1


RA = ReArrangement(FM,m1,m2,n1,n2)

S  = svd(RA)

σ  = S.S[1]
fy = sqrt(σ)*S.U[:,1]
fx = sqrt(σ)*S.Vt[1,:]

Fx = reshape(fx,lx1,lx1)
Fy = reshape(fy,ly1,ly1)

Fld_rec = kron(diag(Fy)',diag(Fx))        # Reconstructed field


close("all")
cm2   = get_cmap("hot");
h1    = figure(num=1,figsize=[16.0,6.0])
axs   = h1.subplots(1,3,sharey=true)
pcm1  = axs[1].pcolormesh(X2D,Y2D,Fld)
pcm1.set_cmap(cm2)

pcm2  = axs[2].pcolormesh(X2D,Y2D,Fld_rec)
pcm2.set_cmap(cm2)

pcm3  = axs[3].pcolormesh(X2D,Y2D,(Fld .- Fld_rec))
pcm3.set_cmap(cm2)

#ax3   = h3.gca()
#ax3.invert_yaxis()
#ax3.set_ylabel("t",fontsize=lafs)
#ax3.set_xlabel("x",fontsize=lafs)
# cb    = colorbar(orientation="vertical")


@printf "Done"
















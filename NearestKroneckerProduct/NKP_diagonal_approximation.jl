#!/bin/julia

println("Nearest Kronecker Product Approximation")

using LinearAlgebra
using PolynomialBases
using PyPlot
using Random
using Printf

include("../TensorOPs/TensorOPs.jl")

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


Nx    = 8
lx1   = Nx+1
Bx    = LobattoLegendre(Nx)
Vx    = ones(Float64,lx1,1)
Nxpf  = 100
XfRef = LinRange(-1.0,1.0,Nxpf)
IntpX = interpolation_matrix(XfRef,Bx.nodes,Bx.baryweights)

Ny    = 12
ly1   = Ny+1
By    = LobattoLegendre(Ny)
Vy    = ones(Float64,ly1,1)
Nypf  = 100
YfRef = LinRange(-1.0,1.0,Nypf)
IntpY = interpolation_matrix(YfRef,By.nodes,By.baryweights)

Nz    = 8
lz1   = Nz+1
Bz    = LobattoLegendre(Nz)
Vz    = ones(Float64,lz1,1)
Nzpf  = 100
ZfRef = LinRange(-1.0,1.0,Nzpf)
IntpZ = interpolation_matrix(ZfRef,Bz.nodes,Bz.baryweights)

xs    = 0.0
xe    = 1.0
ys    = 1.0
ye    = 2.0
zs    = 0.0
ze    = 1.0

Wk2   = zeros(Float64,Nxpf,Nypf)
Wk3   = zeros(Float64,Nxpf,Nypf,Nzpf)
one2  = ones(Float64,1,1)
one3  = ones(Float64,1,1,1)

# X-direction
X     = zeros(Float64,lx1)
for i in 1:lx1
  X[i] = map_from_canonical(Bx.nodes[i],xs,xe,Bx)
end
XM    = reshape(X,lx1,1)
X2D   = zeros(Float64,lx1,ly1)
Xf2D  = zeros(Float64,Nxpf,Nypf)

X3D   = zeros(Float64,lx1,ly1,lz1)
Xf3D  = zeros(Float64,Nxpf,Nypf,Nzpf)

tensorOP2D!(X2D,one2,XM,Vy,Wk2)
tensorOP3D!(X3D,one3,XM,Vy,Vz,Wk3)

tensorOP2D!(Xf2D,X2D,IntpX,IntpY,Wk2)
tensorOP3D!(Xf3D,X3D,IntpX,IntpY,IntpZ,Wk3)

# Y-direction
Y     = zeros(Float64,ly1)
for i in 1:ly1
  Y[i] = map_from_canonical(By.nodes[i],ys,ye,By)
end
#Y2D   = kron(Y',Vx)
#Yf2D  = IntpX*Y2D*(IntpY');

YM    = reshape(Y,ly1,1)
Y2D   = zeros(Float64,lx1,ly1)
Yf2D  = zeros(Float64,Nxpf,Nypf)

Y3D   = zeros(Float64,lx1,ly1,lz1)
Yf3D  = zeros(Float64,Nxpf,Nypf,Nzpf)

tensorOP2D!(Y2D,one2,Vx,YM,Wk2)
tensorOP3D!(Y3D,one3,Vx,YM,Vz,Wk3)

tensorOP2D!(Yf2D,Y2D,IntpX,IntpY,Wk2)
tensorOP3D!(Yf3D,Y3D,IntpX,IntpY,IntpZ,Wk3)

# Z-direction
Z     = zeros(Float64,lz1)
for i in 1:lz1
  Z[i] = map_from_canonical(By.nodes[i],zs,ze,Bz)
end
ZM    = reshape(Z,lz1,1)
#Z2D   = zeros(Float64,lx1,ly1)
#Zf2D  = zeros(Float64,Nxpf,Nypf)

Y3D   = zeros(Float64,lx1,ly1,lz1)
Yf3D  = zeros(Float64,Nxpf,Nypf,Nzpf)

#tensorOP2D!(Y2D,one2,Vx,YM,Wk2)
tensorOP3D!(Y3D,one3,Vx,Vy,ZM,Wk3)

#tensorOP2D!(Yf2D,Y2D,IntpX,IntpY,Wk2)
tensorOP3D!(Yf3D,Y3D,IntpX,IntpY,IntpZ,Wk3)


#rng   = MersenneTwister(1234)
TsFld = zeros(Float64,lx1,ly1)
ErFld = zeros(Float64,lx1,ly1)
Fld   = zeros(Float64,lx1,ly1)
for ci in CartesianIndices(Fld)
  x         = X2D[ci]
  y         = Y2D[ci]
  ErFld[ci] = 1.0*(rand(rng) - 0.5)
  TsFld[ci] = 1.0*cos(8.0*π*x)*sin(5.5*π*y)
  Fld[ci]   = TsFld[ci] + ErFld[ci]
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

TsFldf = IntpX*TsFld*(IntpY')
ErFldf = IntpX*ErFld*(IntpY')
Fldf   = IntpX*Fld*(IntpY')
Fldf_r = IntpX*Fld_rec*(IntpY')
Fmin   = minimum(Fldf[:])
Fmax   = maximum(Fldf[:])


cm2    = get_cmap("PuOr_r"); # RdBu_r
h1,axs = subplots(1,3,sharey=true,figsize=[16.0,6.0],layout="constrained" )

# Underlying Tensor Field
pcm1   = axs[1].pcolormesh(Xf2D,Yf2D,TsFldf,vmin=Fmin,vmax=Fmax)
pcm1.set_cmap(cm2)
cb1    = colorbar(pcm1,location="left")
axs[1].set_title("Tensor Field")

# Total Field
pcm2   = axs[2].pcolormesh(Xf2D,Yf2D,Fldf,vmin=Fmin,vmax=Fmax)
pcm2.set_cmap(cm2)
#cb1    = colorbar(pcm1,location="left")
axs[2].set_title("Total Field")

# Reconstructed Field
pcm3   = axs[3].pcolormesh(Xf2D,Yf2D,Fldf_r,vmin=Fmin,vmax=Fmax)
pcm3.set_cmap(cm2)
#cb2    = colorbar(pcm2,location="top")
axs[3].set_title("Reconstructed Tensor Field")

# pcm4   = axs[4].pcolormesh(Xf2D,Yf2D,(Fldf .- Fldf_r),vmin=Fmin,vmax=Fmax)
# pcm4.set_cmap(cm2)



@printf "Done"
















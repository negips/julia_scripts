# Testing IRAM implementation
println("Testing Implicitly Restarted BiOrthogonal method")

using LinearAlgebra
using Random
#using Pseudospectra
using Printf
using PyPlot

include("BiOrthoUpd.jl")
include("BiOrthoIRst.jl")
include("ExplicitShiftedQR.jl")
include("IRBiOrtho.jl")
include("BulgeChase.jl")
include("HessenbergReduction.jl")

close("all")

rng = MersenneTwister(1234)

vt = ComplexF64
#vt = Float64

zro   = vt(0.0)

n = 20     # Matrix size

λ  = randn(rng,vt,n)
λm = diagm(0 => λ)

U = randn(rng,vt,n,n)
u = qr(U)
UQ = u.Q

A = inv(UQ)*λm*UQ
#A = inv(U)*λm*U

#A = Pseudospectra.grcar(n)
λ = eigvals(A);
AH = A';

λr = real.(λ)
ind = sortperm(λr,rev=true)

Nev   = 5            # Number of eigenvalues to calculate
EKryl = 2            # Additional size of Krylov space
Lk = Nev + EKryl     # Total Size of Krylov space    

Vg    = zeros(vt,n,Lk+1)   # Right Krylov space
Wg    = zeros(vt,n,Lk+1)   # Left  Krylov space
Hv    = zeros(vt,Lk+1,Lk)
Hw    = zeros(vt,Lk+1,Lk)  #

u     = randn(rng,vt,n)
w     = randn(rng,vt,n)
#w     = copy(u)

nk   = 0

ngs = 1
mi        = 1     # Major iteration
nk,mi,ifc = IRBiOrtho!(Vg,Wg,Hv,Hw,u,w,nk,Lk,mi,Nev,ngs)


# Major Iterations
mimax = 2 
while mi < mimax
  global V,W,Hv,Hw,nk,mi

  Av           = A*Vg[:,nk]
  AHw          = AH*Wg[:,nk]

  nk,mi,ifconv = IRBiOrtho!(Vg,Wg,Hv,Hw,Av,AHw,nk,Lk,mi,Nev,ngs)

  if (ifconv)
    break
  end  
 
end  

inp   = Wg'*Vg;

inp_red = inp[1:nk,1:nk];
ortho   = norm(inp_red - I)

println("Subspace Orthogonality: $ortho")

Hu    = copy(Hv[1:Nev,1:Nev])
Hl    = copy(Hw[1:Nev,1:Nev])
λu    = eigvals(Hu)
λl    = eigvals(Hl)
ind2  = sortperm(real.(λu),rev=true)   
ind3  = sortperm(real.(λl),rev=true)   

display([λu[ind2] λ[ind[1:Nev]]])


plot(imag.(λ),real.(λ),linestyle="none", marker="o")
plot(imag(λu),real(λu),linestyle="none", marker="s")


H0 = copy(Hw[1:Nev,1:Nev])
H1 = copy(H0)

#VV = UpperHessenbergReduction(H1)

#k = 1
#v = copy(H0[:,k])
#v[1:k] = zeros(vt,k)
#w = zeros(vt,1,Nev)
#w[k+1] = H0[k,k+1]
#v[k+1] = 1.0/w[k+1]
#β = w*H0[:,k]
#
#Q = I - v*w./β

Vi,Wi = LowerHessenbergtoTriDiagonal!(H1)

I0 = Matrix{vt}(I,Nev,Nev)

v1 = copy(I0)
v1[:,2] = Vi[:,2]


H2 = copy(Hv[1:Nev,1:Nev])
H3 = copy(H2)
V2i,W2i = UpperHessenbergtoTriDiagonal!(H3)


ν = eigvals(H3)
#Hbulge,qbulge = CreateBulge(H3,1,ν,1)
#h1    = copy(Hbulge)
#v1,w1 = LowerHessenbergtoTriDiagonal!(h1)

x = zeros(vt,Nev)

h1 = copy(H3) - ν[1]*I
h1[1,1] = h1[1,1] #- ν[1]
x[1] = zro # H3[1,1] - ν[1]
x[2] = h1[2,1]

y = transpose(h1[1,:])

α = y*x
y = y./α

β = y*h1[:,1] 
Q = I - x*y./β
Qinv = copy(Q)

h2 = Q*h1*Qinv
h3 = copy(h2)

v2,w2 = LowerHessenbergtoTriDiagonal!(h2)


println("Done.")














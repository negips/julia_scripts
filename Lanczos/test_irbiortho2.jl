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

close("all")

rng = MersenneTwister(1234)

vt = ComplexF64
#vt = Float64

zro   = vt(0.0)

n = 100     # Matrix size

λ  = randn(rng,vt,n)
#λ .= λ.^3
λm = diagm(0 => λ)

U = randn(rng,vt,n,n)
u = qr(U)
UQ = u.Q

A = inv(UQ)*λm*UQ
#A = inv(U)*λm*U

#A = Pseudospectra.grcar(n)
AH = A';

λr = real.(λ)
ind = sortperm(λr,rev=true)

Nev   = 5             # Number of eigenvalues to calculate
EKryl = 12             # Additional size of Krylov space
Lk = Nev + EKryl   # Total Size of Krylov space    

V     = zeros(vt,n,Lk+1)   # Right Krylov space
W     = zeros(vt,n,Lk+1)   # Left  Krylov space
Hv    = zeros(vt,Lk+1,Lk)
Hw    = zeros(vt,Lk+1,Lk)
γv    = zeros(vt,Lk+1)
γw    = zeros(vt,Lk+1)

u     = randn(rng,vt,n)
w     = randn(rng,vt,n)

nk   = 0

ifconv = false

BiOrthoUpd!(u,w,V,W,γv,γw,nk)

nk       = nk+1
V[:,nk]  = u
W[:,nk]  = w

# Major Iterations
mi = 1
while mi < 18
  global V,W,Hv,Hw,nk,mi

  Av           = A*V[:,nk]
  AHw          = AH*W[:,nk]

  ngs = 2
  nk,mi = IRBiOrtho!(V,W,Hv,Hw,Av,AHw,nk,Lk,mi,Nev,ngs)
 
end  

Hu    = copy(Hv[1:Nev,1:Nev])
λu    = eigvals(Hu)
ind2  = sortperm(real.(λu),rev=true)   

display([λu[ind2] λ[ind[1:Nev]]])
println("Done.")






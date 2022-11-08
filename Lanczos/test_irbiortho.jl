# Testing IRAM implementation
println("Testing Implicitly Restarted BiOrthogonal method")

using LinearAlgebra
using Random
#using Pseudospectra
using Printf
using PyPlot

include("BiorthoUpd.jl")
#include("ArnIRst.jl")
#include("ExplicitShiftedQR.jl")
#include("IRAM.jl")

close("all")

rng = MersenneTwister(1234)

#vt = ComplexF64
vt = Float64

zro   = vt(0.0)

n = 20     # Matrix size

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
EKryl = 5             # Additional size of Krylov space
LKryl = Nev + EKryl   # Total Size of Krylov space    

Q     = zeros(vt,n,LKryl+1)   # Right Krylov space
P     = zeros(vt,n,LKryl+1)   # Left  Krylov space
Hv    = zeros(vt,LKryl+1,LKryl)
Hw    = zeros(vt,LKryl+1,LKryl)
γv    = zeros(vt,LKryl+1)
γw    = zeros(vt,LKryl+1)

u     = randn(rng,vt,n)
w     = randn(rng,vt,n)

nkryl   = 0

ifconv = false

BiorthoUpd!(u,w,Q,P,γv,γw,nkryl)

nkryl       = nkryl+1
Q[:,nkryl]  = u
P[:,nkryl]  = w

while nkryl<LKryl+1
  global Q,P,Hv,Hw,nkryl

  Au           = A*Q[:,nkryl]
  AHw          = AH*P[:,nkryl]

  BiorthoUpd!(Au,AHw,Q,P,γv,γw,nkryl)

  k = nkryl
  Hv[:,k] = γv
  Hw[:,k] = γw

  Q[:,k+1]  = Au
  P[:,k+1]  = AHw
  nkryl     = nkryl+1

end  
  






# Testing IRAM implementation
println("Testing Implicitly Restarted Lanczos method")

using LinearAlgebra
using Random
#using Pseudospectra
using Printf
using PyPlot

include("LanczosUpd.jl")
#include("ArnIRst.jl")
#include("ExplicitShiftedQR.jl")
#include("IRAM.jl")

close("all")

rng = MersenneTwister(1234)

#vt = ComplexF64
vt = Float64

zro   = vt(0.0)

n = 10     # Matrix size

λ  = randn(rng,vt,n)
#λ .= λ.^3
λm = diagm(0 => λ)

U = randn(rng,vt,n,n)
u = qr(U)
UQ = u.Q

A = inv(UQ)*λm*UQ

#A = Pseudospectra.grcar(n)
AH = A';

λr = real.(λ)
ind = sortperm(λr,rev=true)

Nev   = 4             # Number of eigenvalues to calculate
EKryl = 0            # Additional size of Krylov space
LKryl = Nev + EKryl   # Total Size of Krylov space    

Q     = zeros(vt,n,LKryl+1)   # Left Krylov space
P     = zeros(vt,n,LKryl+1)   # Right Krylov space
Tj    = zeros(vt,LKryl,LKryl)

u     = randn(rng,vt,n)
w     = randn(rng,vt,n)

ngs     = 2       # Number of Gram-Schmidt
nkryl   = 0

α       = zro 
β       = zro 
δ       = zro 

γ       = [α; β; δ]

ifconv = false

LanczosUpd!(u,w,Q,P,γ,nkryl)
display(γ)

nkryl       = nkryl+1
Q[:,nkryl]  = u
P[:,nkryl]  = w
# Hes[1,1] = β
# Hes[1,2] = δ

while nkryl<LKryl+1
  global P,Q,Hes,nkryl
  global γ,u,w
  local α,β,δ

  Au           = A*Q[:,nkryl]
  AHw          = AH*P[:,nkryl]

  LanczosUpd!(Au,AHw,Q,P,γ,nkryl)
  α   = γ[1]
  β   = γ[2]
  δ   = γ[3]

  k = nkryl 
  if k==LKryl
    Tj[k,k]      = α
  else
    Tj[k,k]      = α
    Tj[k+1,k]    = δ
    Tj[k,k+1]    = β
  end  

  nkryl = nkryl+1
  if (nkryl<=LKryl+1)
    Q[:,nkryl]     = Au
    P[:,nkryl]     = AHw
  end 

end  
  






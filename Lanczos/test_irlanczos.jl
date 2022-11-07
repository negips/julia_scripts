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

vt = ComplexF64
#vt = Float64

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
AT = A';

λr = real.(λ)
ind = sortperm(λr,rev=true)

Nev   = 4             # Number of eigenvalues to calculate
EKryl = 0            # Additional size of Krylov space
LKryl = Nev + EKryl   # Total Size of Krylov space    

Q     = zeros(vt,n,LKryl+1)   # Left Krylov space
P     = zeros(vt,n,LKryl+1)   # Right Krylov space
Hes   = zeros(vt,LKryl+1,LKryl)

u     = randn(rng,vt,n)
w     = randn(rng,vt,n)

ngs     = 2       # Number of Gram-Schmidt
nkryl   = 0

α       = zro 
β       = zro 
δ       = zro 

γ       = [α; δ; β]

ifconv = false

LanczosUpd!(Q,P,u,w,γ,nkryl)
display(γ)

nkryl       = nkryl+1
Q[:,nkryl]  = u
P[:,nkryl]  = w
# Hes[1,1] = β
# Hes[1,2] = δ

for i = 1:1 #Nev
  global P,Q,Hes,nkryl
  global γ,u,w
  local α,β,δ

  u           = A*u
  w           = AT*w

  LanczosUpd!(Q,P,u,w,γ,nkryl)
  nkryl            = nkryl+1
  α   = γ[1]
  δ   = γ[2]
  β   = γ[3]

  k = nkryl -1 
  if k==1
    Hes[k,k]      = δ
    Hes[k+1,k]    = β
    Hes[k+1,k]    = β

  elseif k==LKryl
    Hes[k-1,k]    = α
    Hes[k,k]      = δ
  else
    Hes[k-1,k]    = α
    Hes[k,k]      = δ
    Hes[k+1,k]    = β
  end  

  if (nkryl<LKryl)
    Q[:,nkryl]     = u
    P[:,nkryl]     = w
  end 

   

end  
  






# Testing IRAM implementation
println("Testing Implicitly Restarted Lanczos method")

using LinearAlgebra
using Random
using Pseudospectra
using Printf
using PyPlot

include("ArnUpd.jl")
include("ArnIRst.jl")
include("ExplicitShiftedQR.jl")
include("IRAM.jl")

close("all")

rng = MersenneTwister(1234)

vt = ComplexF64
#vt = Float64

n = 100     # Matrix size

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

Nev   = 20             # Number of eigenvalues to calculate
EKryl = 36            # Additional size of Krylov space
LKryl = Nev + EKryl   # Total Size of Krylov space    

Q     = zeros(vt,n,LKryl+1)   # Left Krylov space
P     = zeros(vt,n,LKryl+1)   # Right Krylov space
Hes   = zeros(vt,LKryl+1,LKryl)

u     = randn(rng,vt,n)
w     = randn(rng,vt,n)

ngs     = 2       # Number of Gram-Schmidt
nkryl   = 0
β       = norm(u)      # ||u0||
ϵ       = v'*w/β0
α       = 0.0

ifconv = false

for i = 1:Nev
  global P,Q,Hes,nkryl
  local h,α,β,ϵ,w,u

  q  = u./β
  p  = w./ϵ

  Q[:,i] = q
  P[:,i] = p
 
    u  = A*q
    w  = AT*p
  
  if (i>1)
    u = u .- Q[:,i-1]*ϵ
    w = w .- P[:,i-1]*β
  end

  α = p'*u
  u = u .- q*α
  w = w .- p*α

  β = norm(u)
  ϵ = u'*w./β


end  
  


#  α,β,ϵ  = LanczosUpd(W,V,y,x,nkryl,ngs)
#  Hes[1:nkryl,nkryl] = h
#  Hes[nkryl+1,nkryl] = β
#  V[:,nkryl+1]       = r
#  nkryl              = nkryl + 1


#while ~ifconv

## Major Iterations
#for mi in 1:500
#  global V,Hes,nkryl,ifconv
#  local U,G
#  local β
#
#  for i in Nev+1:LKryl
#
#    local h,β,r,v
#
#    r = V[:,nkryl]
#    v = A*r
#
#    V,Hes,nkryl,β,mi2 = IRAM!(V,Hes,Bg,v,nkryl,LKryl,mi,Nev,ngs)
#
#  end
#
#  β   = abs(Hes[Nev+1,Nev])
#
#  if β < 1.0e-12
#    break
#  end
# 
#end  

#Ht = Hes[1:Nev,1:Nev]
#F = eigen(Ht)
#ind2 = sortperm(real.(F.values),rev=true)

#display(F.values[ind2])

#plot(real.(F.values),imag.(F.values),linestyle="none",marker="o")




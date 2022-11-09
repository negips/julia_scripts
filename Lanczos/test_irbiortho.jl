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

V     = zeros(vt,n,LKryl+1)   # Right Krylov space
W     = zeros(vt,n,LKryl+1)   # Left  Krylov space
Hv    = zeros(vt,LKryl+1,LKryl)
Hw    = zeros(vt,LKryl+1,LKryl)
γv    = zeros(vt,LKryl+1)
γw    = zeros(vt,LKryl+1)

u     = randn(rng,vt,n)
w     = randn(rng,vt,n)

nkryl   = 0

ifconv = false

BiOrthoUpd!(u,w,V,W,γv,γw,nkryl)

nkryl       = nkryl+1
V[:,nkryl]  = u
W[:,nkryl]  = w
#
#while nkryl<Nev+1
#  global V,W,Hv,Hw,nkryl
#
#  Av           = A*V[:,nkryl]
#  AHw          = AH*W[:,nkryl]
#
#  BiOrthoUpd!(Av,AHw,V,W,γv,γw,nkryl)
#
#  k = nkryl
#  Hv[:,k] = γv
#  Hw[:,k] = γw'
#
#  V[:,k+1]  = Av
#  W[:,k+1]  = AHw
#  nkryl     = nkryl+1
#
#end  
  



# Major Iterations
mi = 1
while mi < 2
  global V,W,Hv,Hw,nkryl,mi

  Av           = A*V[:,nkryl]
  AHw          = AH*W[:,nkryl]

  ngs = 2
  nkryl,mi = IRBiOrtho!(V,W,Hv,Hw,Av,AHw,nkryl,LKryl,mi,Nev,ngs)
  println("nkryl=$nkryl, Major Iteration=$mi")
 
end  

#H     = copy(Hv[1:LKryl,1:LKryl])
#μ,nμ  = ArnGetLowerShifts(H,EKryl)
#
#ngs   = 2
#Hs,Q  = ExplicitShiftedQR(H,μ,nμ,ngs)

println("Done.")






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

n = 300     # Matrix size

λ  = randn(rng,vt,n)
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

Nev   = 5            # Number of eigenvalues to calculate
EKryl = 10            # Additional size of Krylov space
Lk = Nev + EKryl   # Total Size of Krylov space    

V     = zeros(vt,n,Lk+1)   # Right Krylov space
W     = zeros(vt,n,Lk+1)   # Left  Krylov space
Hv    = zeros(vt,Lk+1,Lk)
Hw    = zeros(vt,Lk+1,Lk)  #
γv    = zeros(vt,Lk+1)
γw    = zeros(vt,Lk+1)

u     = randn(rng,vt,n)
#w     = randn(rng,vt,n)
w     = copy(u)

nk   = 0

BiOrthoUpd!(u,w,V,W,γv,γw,nk)

nk       = nk+1
V[:,nk]  = u
W[:,nk]  = w

# Major Iterations
mi = 1
mimax = 50 
while mi < mimax
  global V,W,Hv,Hw,nk,mi

  Av           = A*V[:,nk]
  AHw          = AH*W[:,nk]

  ngs = 2
  nk,mi,ifconv = IRBiOrtho!(V,W,Hv,Hw,Av,AHw,nk,Lk,mi,Nev,ngs)

  if (ifconv)
    break
  end  
 
end  

inp   = W'*V;

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

#ekryl = Lk - Nev;
#
#
#H = copy(Hv[1:Lk,1:Lk])
#
##Perform QR operations
#ngs     = 2
#μ,nμ    = GetLowerShifts(H,ekryl)
#
#Hsv,Qv  = ExplicitShiftedQR(H,μ,nμ,ngs)   # Hsv = Q'*H*Q; H = Q*Hsv*Q'
#
#Hsw     = Qv'*Hw[1:Lk,1:Lk]*Qv 
#
#Vp      = V[:,1:Lk]*Qv
#Wp      = W[:,1:Lk]*Qv
#
#vk1     = Hv[Lk+1,Lk]*V[:,Lk+1];
#wk1     = Hw[Lk+1,Lk]*W[:,Lk+1];
#
#ek      = zeros(vt,Lk)
#ek[Lk]  = vt(1.0)
#
#ektQ    = transpose(Qv[Lk,:])
#
#rMv     = vk1*ektQ      # residual matrix v
#rMw     = wk1*ektQ      # residual matrix w
#
#rv2     = rMv[:,Nev]    # residual vector
#rw2     = rMw[:,Nev]    # residual vector


println("Done.")










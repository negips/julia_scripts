# Testing IRAM implementation
println("Testing Implicitly Restarted BiOrthogonal method")

using LinearAlgebra
using Random
using Pseudospectra
using Printf
using PyPlot

include("BiOrthoUpd.jl")
include("BiOrthoIRst.jl")
include("ExplicitShiftedQR.jl")
include("IRBiOrtho.jl")
include("BulgeChase.jl")
include("HessenbergReduction.jl")
include("NegiEig.jl")
include("ImplicitLR.jl")
include("SetZero.jl")

close("all")

rng = MersenneTwister(1234)

vt = ComplexF64
#vt = Float64

zro   = vt(0.0)
one   = vt(1.0)

n = 40     # Matrix size

λ  = randn(rng,vt,n)
λm = diagm(0 => λ)

U = randn(rng,vt,n,n)
u = qr(U)
UQ = u.Q

A = inv(UQ)*λm*UQ
#A = inv(U)*λm*U

A = Pseudospectra.grcar(n)
λ = eigvals(A);
AH = A';

λr = real.(λ)
ind = sortperm(λr,rev=true)

Nev   = 10            # Number of eigenvalues to calculate
EKryl = 10            # Additional size of Krylov space
Lk = Nev + EKryl     # Total Size of Krylov space    

Vg    = zeros(vt,n,Lk+1)   # Right Krylov space
Wg    = zeros(vt,n,Lk+1)   # Left  Krylov space
Hv    = zeros(vt,Lk+1,Lk)
Hw    = zeros(vt,Lk+1,Lk)  #

#v     = randn(rng,vt,n)
#w     = randn(rng,vt,n)

v     = randn(vt,n)
w     = randn(vt,n)

nk   = 0

ngs = 2
mi        = 1     # Major iteration
nk,mi,ifc = IRBiOrtho2!(Vg,Wg,Hv,Hw,v,w,nk,Lk,mi,Nev,ngs)

# Major Iterations
mimax = 65 
while mi < mimax
  global Vg,Wg,Hv,Hw,nk,mi

  miold        = mi
  Av           = A*Vg[:,nk]
  AHw          = AH*Wg[:,nk]

  nk,mi,ifconv = IRBiOrtho2!(Vg,Wg,Hv,Hw,Av,AHw,nk,Lk,mi,Nev,ngs)

  if (mi>miold)

#    BiorthoReortho!(Vg,Wg,nk,2)
#    BiorthoReortho2!(Vg,Wg,nk,2)
#    BiorthoReortho3!(Vg,Wg,nk,2)

    inp   = Wg'*Vg;
    inp_red = inp[1:nk,1:nk];
    ortho   = norm(inp_red - I)
    
    println("Subspace BiOrthogonality: $mi: $ortho")

    epsi  = inp - I;
#    Ipert = (I - epsi)';            # Using this as approximate inverse
#    Wg    = Wg*Ipert

#    Wg[:,1:nk] = Wg[:,1:nk]*(inv(inp_red)')
#
#    inp   = Wg'*Vg;
#    inp_red = inp[1:nk,1:nk];
#    ortho   = norm(inp_red - I)
#   
#    println("Subspace BiOrthogonality (corrected): $mi: $ortho")

#    Vg .= Vg .- Vg*(inp - I)
#    Wg .= Wg .- Wg*(inp' - I)

  end  


  if (ifconv)
    break
  end  
 
end  

H = Hv[1:Nev,1:Nev]
θ = eigvals(H)

inp   = Wg'*Vg;

inp_red = inp[1:nk,1:nk];
ortho   = norm(inp_red - I)

println("Subspace BiOrthogonality1: $ortho")

close("all")

#plot(imag.(λ),real.(λ),linestyle="none",marker="o")
#plot(imag.(θ),real.(θ),linestyle="none",marker="*",markersize=8)

plot(real.(λ),imag.(λ),linestyle="none",marker="o")
plot(real.(θ),imag.(θ),linestyle="none",marker="*",markersize=8)

#BiorthoReortho2!(Vg,Wg,nk,2)

#Verr = Vg*(inp - I)
#Vg2  = Vg .- Verr
inp2 = Wg'*Vg;
inp_red2 = inp2[1:nk,1:nk]

ortho2   = norm(inp_red2 - I)
println("Subspace BiOrthogonality2: $ortho2")

BiorthoReortho2!(Vg,Wg,nk,2)

#Verr2 = Vg2*(inp2 - I)
#Vg3   = Vg2 .- Verr2
inp3  = Wg'*Vg;

inp_red3 = inp3[1:nk,1:nk]
ortho3   = norm(inp_red3 - I)
println("Subspace BiOrthogonality3: $ortho3")



println("Done.")














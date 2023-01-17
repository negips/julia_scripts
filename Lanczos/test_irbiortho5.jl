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
include("NegiEig.jl")
include("ImplicitLR.jl")
include("SetZero.jl")

close("all")

rng = MersenneTwister(1234)

vt = ComplexF64
#vt = Float64

zro   = vt(0.0)
one   = vt(1.0)

n = 15     # Matrix size

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

Nev   = 3            # Number of eigenvalues to calculate
EKryl = 5            # Additional size of Krylov space
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

H = Hv[1:Lk,1:Lk]
SetZero!(H,1.0e-13)
Hsv = copy(H)
μ,nμ = GetLowerShifts(Hsv,EKryl);
vi,wi = ImplicitLRSeq!(Hsv,μ,nμ);   # wi*H*vi = Hsv
SetZero!(Hsv,1.0e-12)


Vp   = Vg[:,1:Lk]*vi;
Wp   = Wg[:,1:Lk]*(wi');


av = Vg[:,Lk+1]*Hv[Lk+1,Lk]
aw = Wg[:,Lk+1]*Hw[Lk+1,Lk]

ekt = zeros(vt,1,Lk)
ekt[Lk] = one

rvM = av*ekt*vi
rwM = aw*ekt*(wi')

vv  = Hsv[Nev+1,Nev]*Vp[:,Nev] .+ rvM[:,Nev] 
ww  = Hsv'[Nev+1,Nev]*Wp[:,Nev] .+ rwM[:,Nev] 


H2 = copy(Hsv');
v2,w2 = ImplicitLRSeq!(H2,adjoint.(μ),nμ);   # wi*H*vi = Hsv

tst1 = vi*(w2')

tst2 = wi'*v2



inp   = Wg'*Vg;

inp_red = inp[1:nk,1:nk];
ortho   = norm(inp_red - I)

println("Subspace Orthogonality: $ortho")


println("Done.")














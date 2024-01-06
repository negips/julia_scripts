# Testing IRAM implementation
println("Testing Implicitly Restarted BiOrthogonal method")

using LinearAlgebra
using Random
using Pseudospectra
using Printf
using PyPlot,PyCall

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

ifsavefig = false

rng = MersenneTwister(1236)   # paper: 1236

vt = ComplexF64
#vt = Float64

zro   = vt(0.0)
one   = vt(1.0)

lafs  = 16


n = 50     # Matrix size

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

λi  = imag.(λ)
indi = sortperm(λi,rev=true)


Nev   = 10            # Number of eigenvalues to calculate
EKryl = 10            # Additional size of Krylov space
Lk    = Nev + EKryl   # Total Size of Krylov space    

λ2    = λ[indi[1:Nev]]

Vg    = zeros(vt,n,Lk+1)   # Right Krylov space
Wg    = zeros(vt,n,Lk+1)   # Left  Krylov space
Hv    = zeros(vt,Lk+1,Lk)
Hw    = zeros(vt,Lk+1,Lk)  #

v     = randn(rng,vt,n)
w     = randn(rng,vt,n)

#v     = randn(vt,n)
#w     = randn(vt,n)

nk   = 0

ngs = 2
mi        = 1     # Major iteration
nk,mi,ifc = IRBiOrtho2!(Vg,Wg,Hv,Hw,v,w,nk,Lk,mi,Nev,ngs)



# Major Iterations
mimax = 25

λerr  = zeros(vt,mimax,Nev)

while mi < mimax
  global Vg,Wg,Hv,Hw,nk,mi

  miold        = mi
  Av           = A*Vg[:,nk]
  AHw          = AH*Wg[:,nk]

  nk,mi,ifconv = IRBiOrtho2!(Vg,Wg,Hv,Hw,Av,AHw,nk,Lk,mi,Nev,ngs)

  if (mi>miold)

    inp   = Wg'*Vg;
    inp_red = inp[1:nk,1:nk];
    ortho   = norm(inp_red - I)
    
    println("Subspace BiOrthogonality: $mi: $ortho")

    epsi  = inp - I;
    Ipert = (I - epsi)';            # Using this as approximate inverse
    Wg    = Wg*Ipert

#    Wg[:,1:nk] = Wg[:,1:nk]*(inv(inp_red)')

    inp   = Wg'*Vg;
    inp_red = inp[1:nk,1:nk];
    ortho   = norm(inp_red - I)
#   
    println("Subspace BiOrthogonality (corrected): $mi: $ortho")
#

    Tred = Hv[1:Nev,1:Nev]
    evs  = eigvals(Tred)
    ϕi   = imag.(evs)
    indi2= sortperm(ϕi,rev=true)
      
    ϕ2   = evs[indi2]

    eig_err = λ2 .- ϕ2
    λerr[mi-1,:]  = transpose(eig_err);

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

#matplotlib.rcParams["xtick.labelsize"] = 16
#PyCall.PyDict(matplotlib["rcParams"])["markers.fillstyle"] = "full"
#PyCall.PyDict(matplotlib["rcParams"])["axes.labelsize"] = 16;
#PyCall.PyDict(matplotlib["rcParams"])["axes.grid"] = true;
#PyCall.PyDict(matplotlib["rcParams"])["grid.linestyle"] = ":";

PyCall.PyDict(matplotlib["rcParams"])["xtick.labelsize"] = 12;
PyCall.PyDict(matplotlib["rcParams"])["ytick.labelsize"] = 12;

h1  = figure(num=1,figsize=[8.,8.]);
semilogy(abs.(λerr[1:mi-1,:]))
ax1 = h1.axes[1]
ax1.set_ylabel(L"||λ-λ_{L}||", fontsize=lafs)
ax1.set_xlabel("i", fontsize=lafs)

if (ifsavefig)
  savefig("eigen_error.eps")
end  

#plot(imag.(λ),real.(λ),linestyle="none",marker="o")
#plot(imag.(θ),real.(θ),linestyle="none",marker="*",markersize=8)

PyCall.PyDict(matplotlib["rcParams"])["xtick.labelsize"] = 12;
PyCall.PyDict(matplotlib["rcParams"])["ytick.labelsize"] = 12;

h2  = figure(num=2,figsize=[8.,8.]);
plot(real.(λ),imag.(λ),linestyle="none",marker="o",label="Grcar Spectra")
plot(real.(θ),imag.(θ),linestyle="none",marker="s",markersize=8,markeredgewidth=2,fillstyle="none",label="Restarted Lanczos")
ax2 = h2.axes[1]
ax2.set_xlabel(L"λ_{r}", fontsize=lafs)
ax2.set_ylabel(L"λ_{i}", fontsize=lafs)
ax2.legend(fontsize=12)

if (ifsavefig)
  savefig("grcar_eigs.eps")
end  


println("Done.")














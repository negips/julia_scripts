# Testing IRAM implementation
println("Testing Implicitly Restarted BiOrthogonal method")

using LinearAlgebra
using Random
using Pseudospectra
using Printf
using PyPlot,PyCall

include("Module_SchurLanczos/SchurLanczos.jl")
include("SetZero.jl")

close("all")

ifsavefig = false

rng = MersenneTwister(1236)   # paper: 1236

const vt = ComplexF64
const BiOtol = 1.0e-08        # small biorthogonality tolerance
const SStol  = 1.0e-08        # invariant subspace tolerance
const zrotol = 1.0e-12        # Zero element tolerance

#vt = Float64

zro   = vt(0.0)
one   = vt(1.0)

lafs  = 16

n  = 10     # Matrix size

A  = Pseudospectra.grcar(n)
λ  = eigvals(A);
AH = A';

λr  = real.(λ)
ind = sortperm(λr,rev=true)

λi   = imag.(λ)
indi = sortperm(λi,rev=true)


Nev   = 2             # Number of eigenvalues to calculate
EKryl = 3             # Additional size of Krylov space
Lk    = Nev + EKryl   # Total Size of Krylov space
Lk1   = Lk + 1

λ2    = λ[indi[1:Nev]]

Qr    = zeros(vt,n,Lk)   # Right Krylov space
Ql    = zeros(vt,n,Lk)   # Left  Krylov space
Rr    = zeros(vt,Lk,Lk)
Rl    = zeros(vt,Lk,Lk)  
H     = zeros(vt,Lk,Lk)
γl    = zeros(vt,Lk)
γr    = zeros(vt,Lk)

v     = randn(rng,vt,n)
w     = randn(rng,vt,n)
x     = zeros(vt,n)     # Work vector

nk    = 0

ngs   = 2
mi    = 1     # Major iteration
#nk,mi,ifc = IRBiOrtho2!(Vg,Wg,Hv,Hw,v,w,nk,Lk,mi,Nev,ngs)

k         = 0
Qlv,Rlv   = SchurLanczos.GetQRviews(Ql,Rl,k)
Qrv,Rrv   = SchurLanczos.GetQRviews(Qr,Rr,k)
γlv       = view(γl,1:k)
γrv       = view(γr,1:k)

αl,αr = SchurLanczos.BiOrthoScale!(w,v,x)
SchurLanczos.UpdateQR!(Ql,Rl,w,k,ngs)
SchurLanczos.UpdateQR!(Qr,Rr,v,k,ngs)

include("OneStep.jl")

#k        += 1
#Qlv,Rlv   = SchurLanczos.GetQRviews(Ql,Rl,k)
#Qrv,Rrv   = SchurLanczos.GetQRviews(Qr,Rr,k)
#γlv       = view(γl,1:k)
#γrv       = view(γr,1:k)
#
#w1        = copy(w)
#v1        = copy(v)
#
#v         = A*v
#w         = A'*w
#
#SchurLanczos.BiOrthogonalize!(γlv,γrv,w,v,Qlv,Qrv,Rlv,Rrv,ngs)
#αl,αr = SchurLanczos.BiOrthoScale!(w,v)
#SchurLanczos.UpdateQR!(Ql,Rl,w,k,ngs)
#SchurLanczos.UpdateQR!(Qr,Rr,v,k,ngs)



## # Major Iterations
## mimax = 25
## 
## λerr  = zeros(vt,mimax,Nev)
## 
## while mi < mimax
##   global Vg,Wg,Hv,Hw,nk,mi
## 
##   miold        = mi
##   Av           = A*Vg[:,nk]
##   AHw          = AH*Wg[:,nk]
## 
##   nk,mi,ifconv = IRBiOrtho2!(Vg,Wg,Hv,Hw,Av,AHw,nk,Lk,mi,Nev,ngs)
## 
##   if (mi>miold)
## 
##     inp   = Wg'*Vg;
##     inp_red = inp[1:nk,1:nk];
##     ortho   = norm(inp_red - I)
##     
##     println("Subspace BiOrthogonality: $mi: $ortho")
## 
##     epsi  = inp - I;
##     Ipert = (I - epsi)';            # Using this as approximate inverse
##     Wg    = Wg*Ipert
## 
## #    Wg[:,1:nk] = Wg[:,1:nk]*(inv(inp_red)')
## 
##     inp   = Wg'*Vg;
##     inp_red = inp[1:nk,1:nk];
##     ortho   = norm(inp_red - I)
## #   
##     println("Subspace BiOrthogonality (corrected): $mi: $ortho")
## #
## 
##     Tred = Hv[1:Nev,1:Nev]
##     evs  = eigvals(Tred)
##     ϕi   = imag.(evs)
##     indi2= sortperm(ϕi,rev=true)
##       
##     ϕ2   = evs[indi2]
## 
##     eig_err = λ2 .- ϕ2
##     λerr[mi-1,:]  = transpose(eig_err);
## 
##   end  
## 
## 
##   if (ifconv)
##     break
##   end  
##  
## end  
## 
## H = Hv[1:Nev,1:Nev]
## θ = eigvals(H)
## 
## inp   = Wg'*Vg;
## 
## inp_red = inp[1:nk,1:nk];
## ortho   = norm(inp_red - I)
## 
## println("Subspace BiOrthogonality1: $ortho")
## 
## #matplotlib.rcParams["xtick.labelsize"] = 16
## #PyCall.PyDict(matplotlib["rcParams"])["markers.fillstyle"] = "full"
## #PyCall.PyDict(matplotlib["rcParams"])["axes.labelsize"] = 16;
## #PyCall.PyDict(matplotlib["rcParams"])["axes.grid"] = true;
## #PyCall.PyDict(matplotlib["rcParams"])["grid.linestyle"] = ":";
## 
## PyCall.PyDict(matplotlib["rcParams"])["xtick.labelsize"] = 12;
## PyCall.PyDict(matplotlib["rcParams"])["ytick.labelsize"] = 12;
## 
## h1  = figure(num=1,figsize=[8.,8.]);
## semilogy(abs.(λerr[1:mi-1,:]))
## ax1 = h1.axes[1]
## ax1.set_ylabel(L"||λ-λ_{L}||", fontsize=lafs)
## ax1.set_xlabel("i", fontsize=lafs)
## 
## if (ifsavefig)
##   savefig("eigen_error.eps")
## end  
## 
## #plot(imag.(λ),real.(λ),linestyle="none",marker="o")
## #plot(imag.(θ),real.(θ),linestyle="none",marker="*",markersize=8)
## 
## PyCall.PyDict(matplotlib["rcParams"])["xtick.labelsize"] = 12;
## PyCall.PyDict(matplotlib["rcParams"])["ytick.labelsize"] = 12;
## 
## h2  = figure(num=2,figsize=[8.,8.]);
## plot(real.(λ),imag.(λ),linestyle="none",marker="o",label="Grcar Spectra")
## plot(real.(θ),imag.(θ),linestyle="none",marker="s",markersize=8,markeredgewidth=2,fillstyle="none",label="Restarted Lanczos")
## ax2 = h2.axes[1]
## ax2.set_xlabel(L"λ_{r}", fontsize=lafs)
## ax2.set_ylabel(L"λ_{i}", fontsize=lafs)
## ax2.legend(fontsize=12)
## 
## if (ifsavefig)
##   savefig("grcar_eigs.eps")
## end  


println("Done.")














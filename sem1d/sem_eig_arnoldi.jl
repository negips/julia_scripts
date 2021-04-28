println("Main interface for 1D SEM")

using PolynomialBases
using LinearAlgebra
# using UnicodePlots
#using Plots
using PyPlot,Colors,PyCall
using Arpack

# Include the function files
include("sem_geom.jl")

include("sem_init_ref.jl")

include("AssembleMatrix.jl")

include("AssembleMatrixLesshafft.jl")


close("all")

c0 = 1.0e-10;

L,M,OP,Binv  = AssembleMatrixLesshafft(c0,Geom.cnv,Geom.wlp,Geom.xm1,Geom.bm1,Basis,lx1,nel);

#Labs = abs.(L)

#expL = exp(L)
#eLabs = abs.(expL)

#h1 = figure(num=1,figsize=[8.,6.]);

r,c = size(L)
r1 = 2; r2 = r-1
c1 = 2; c2 = c-1
Lt = L[r1:r2,c1:c2]
Mt = M[r1:r2,c1:c2]
Lnew = inv(Mt)*Lt

#spy(Labs)

Ω = eigvals(Lnew)
# F = eigen(inv(Mt)*Lt)
# Ω = F.values
Lesshafft_ω = 1.0*im*Ω
ωi = imag(Lesshafft_ω)
ωr = real(Lesshafft_ω)

h2 = figure(num=2,figsize=[8.,6.]);
ax1 = gca()
p1 = ax1.plot(ωr,ωi, linestyle="none",marker=".")


# eigs is from the Arpack package
λ, ϕ = eigs(Lnew, nev=16, ncv=50, which=:LR, maxiter=500)

Lesshafft_λ = 1.0*im*λ

λi = imag(Lesshafft_λ)
λr = real(Lesshafft_λ)

p1 = ax1.plot(λr,λi, linestyle="none",marker=".")

ax1.set_xlim(-5.0, 5.0)
ax1.set_ylim(-5.0, 1.0)
grid(true)

println("Done")











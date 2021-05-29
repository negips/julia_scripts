println("Main interface for 1D SEM")

using PolynomialBases
using LinearAlgebra
# using UnicodePlots
#using Plots
using PyPlot,Colors,PyCall

# Include the function files
include("sem_geom.jl")

include("sem_init_ref.jl")

include("AssembleMatrix.jl")

include("AssembleMatrixLesshafft.jl")


close("all")

c0 = 0.0e-10;
ν  = 0.1;   # nu


#L  = AssembleMatrix(c0,Geom.cnv,ν,Geom.wlp,lx1,nel);

L,M,OP  = AssembleMatrixLesshafft(c0,Geom.cnv,Geom.wlp,Geom.xm1,Geom.bm1,Basis,lx1,nel);

Labs = abs.(L)

#expL = exp(L)
#eLabs = abs.(expL)

h1 = figure(num=1,figsize=[8.,6.]);

r,c = size(L)
r1 = 2; r2 = r
c1 = 2; c2 = c
Lt = L[r1:r2,c1:c2]
Mt = M[r1:r2,c1:c2]
spy(Labs)

Ω = eigvals(inv(Mt)*Lt)
# F = eigen(inv(Mt)*Lt)
# Ω = F.values

Lesshafft_ω = 1.0*im*Ω

ωi = imag(Lesshafft_ω)
ωr = real(Lesshafft_ω)

h2 = figure(num=2,figsize=[8.,6.]);
ax1 = gca()
p1 = ax1.plot(ωr,ωi, linestyle="none",marker=".")


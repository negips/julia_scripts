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

include("Sem_QQT.jl")

close("all")

ndof, glnum = Sem_Global_Num(Geom.xm1)

Q,QT  = Sem_QQT(glnum)
v     = sum(Q*QT,dims=2)
#vmult = reshape(v,lx1,nel)
vimult = 1.0 ./v

c0 = 0.0e-10;

L,M,OP,Binv  = AssembleMatrixLesshafft(c0,Geom.cnv,Geom.wlp,Geom.xm1,Geom.bm1,Basis,lx1,nel);


xglob = vimult.*Q*QT*Geom.xm1[:]

Md = diag(M)
Mdinv = 1.0./Md

Minv = diagm(Mdinv)
#Labs = abs.(L)

#expL = exp(L)
#eLabs = abs.(expL)

#h1 = figure(num=1,figsize=[8.,6.]);

r,c = size(L)
r1 = 2; r2 = r
c1 = 2; c2 = c
Lt = L[r1:r2,c1:c2]
Mt = M[r1:r2,c1:c2]
Minvt = Minv[r1:r2,c1:c2]

Lnew = similar(Lt)

mul!(Lnew,Minvt,Lt)
#Lnew = Lt

#spy(Labs)

# Ω = eigvals(Lnew)
#F = eigen(Lnew)
#Ω = F.values
#Lesshafft_ω = 1.0*im*Ω
#ωi = imag(Lesshafft_ω)
#ωr = real(Lesshafft_ω)

h1 = figure(num=1,figsize=[8.,6.]);
ax1 = gca()
#p0 = ax1.plot(ωr,ωi, linestyle="none",marker="o")


# eigs is from the Arpack package
println("Starting IRAM")
λ, ϕ = eigs(Lnew, nev=40, ncv=150, which=:LR, maxiter=500)

Lesshafft_λ = 1.0*im*λ

λi = imag(Lesshafft_λ)
λr = real(Lesshafft_λ)

p1 = ax1.plot(λr,λi, linestyle="none",marker=".",markersize=16)

ax1.set_xlim(-6.0, 6.0)
ax1.set_ylim(-10.0, 1.0)
grid(true)

h2 = figure(num=2,figsize=[8.,6.]);
ax2 = gca()
for i in 1:1
  p2 = ax2.plot(xglob[r1:r2],real.(ϕ[:,i]), linestyle="-")
#  p3 = ax2.plot(xglob[r1:r2],imag.(F.vectors[:,r2-i]), linestyle="-")
end  


println("Done")










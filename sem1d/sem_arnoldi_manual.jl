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
vmult = sum(Q*QT,dims=2)
#vmult = reshape(v,lx1,nel)
vimult = 1.0 ./vmult

c0 = 0.0e-10;

L,M,OP,Binv = AssembleMatrixLesshafft(c0,Geom.cnv,Geom.wlp,Geom.xm1,Geom.bm1,Basis,lx1,nel);

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


#h1 = figure(num=1,figsize=[8.,6.]);
#ax1 = gca()

# eigs is from the Arpack package
println("Starting IRAM")

Nev   = 5
Nkryl = 10
ngs   = 1
B   = Geom.bm1[:]       # Mass matrix

vlen = lx1*nel
V   = zeros(Complex,vlen,Nkryl+1)
H   = zeros(Complex,Nkryl+1,Nkryl)
v0  = rand(lx1*nel,1) + im*rand(lx1*nel,1)
v0  = vimult.*Q*QT*v0               # Make v0 continuous across elements
β   = sqrt(v0'*(B.*v0))             # B-orthogonal
v   = v0/β

V[:,1] = v

# Populate Nev Krylov vectors for the first time
for i in 1:Nev
  global v
  local β
  global V, H
  for e in 1:nel
    j1 = (e-1)*lx1 + 1
    j2 = e*lx1
    v[j1:j2] = OP[:,:,e]*v[j1:j2]
  end
  v  = Q*QT*v
  g  = zeros(Complex,Nkryl+1)
  for j in 1:ngs
    local h
    h  = V[:,1:i]'*v
    v  = v - V[:,1:i]*h
    g[1:i] = g[1:i] .+ h
  end
  vnorm    = sqrt(v'*(B.*v))
  β        = vnorm[1]
  g[i+1]   = β
  H[:,i]   = g
  v        = v/β
  V[:,i+1] = v
end  

  
# (Nkryl-Nev) Krylov vectors
for i in Nev+1:Nkryl
  global v
  local β
  global V, H
  for e in 1:nel
    j1 = (e-1)*lx1 + 1
    j2 = e*lx1
    v[j1:j2] = OP[:,:,e]*v[j1:j2]
  end
  v  = Q*QT*v
  g  = zeros(Complex,Nkryl+1)
  for j in 1:ngs
    local h
    h  = V[:,1:i]'*v
    v  = v - V[:,1:i]*h
    g[1:i] = g[1:i] .+ h
  end
  vnorm    = sqrt(v'*(B.*v))
  β        = vnorm[1]
  g[i+1]   = β
  H[:,i]   = g
  v        = v/β
  V[:,i+1] = v
end  



# λ, ϕ = eigs(Lnew, nev=40, ncv=150, which=:LR, maxiter=500)
# 
# Lesshafft_λ = 1.0*im*λ

# 

# λi = imag(Lesshafft_λ)
# λr = real(Lesshafft_λ)
# 
# p1 = ax1.plot(λr,λi, linestyle="none",marker=".",markersize=16)
# 
# ax1.set_xlim(-6.0, 6.0)
# ax1.set_ylim(-10.0, 1.0)
# grid(true)
# 
# h2 = figure(num=2,figsize=[8.,6.]);
# ax2 = gca()
# for i in 1:1
#   p2 = ax2.plot(xglob[r1:r2],real.(ϕ[:,i]), linestyle="-")
# #  p3 = ax2.plot(xglob[r1:r2],imag.(F.vectors[:,r2-i]), linestyle="-")
# end  
# 
# 
# println("Done")










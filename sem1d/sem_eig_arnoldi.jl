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

manualarnoldi = true

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
Nev = 40
Ncv = 150

if (manualarnoldi)
  vol = sum(Geom.bm1[:])                                    # Volume
  V   = zeros(Complex,lx1*nel,Nev+1)
  H   = zeros(Complex,Nev+1,Nev)
  Hb  = zeros(Complex,Nev,Nev)

  B   = Geom.bm1[:]/vol       # B-orthogonal Weight

  v   = rand(Float64,lx1*nel,1) .+ im.*rand(Float64,lx1*nel,1)    # Starting Vector
  v[1,1] = 0. + im*0.
  vnorm = sqrt(v'*(B.*v))                               # Normalize
  v   = v./vnorm
  V[:,1] = v
 
  v   = reshape(v,lx1,nel)
  for i in 1:Nev
    global v
    local vnorm
#   Apply operator 
    for e in 1:nel
      v[:,e] = OP[:,:,e]*v[:,e];
    end
    v[1,1] = 0. + im*0.
    vg       = Q*QT*v[:]
    vb       = B.*vg
    hij      = V[:,1:i]'*vb         # projections
    H[1:i,i] = hij 
    vp       = V[:,1:i]*hij         # Projected vector      
    vg       = vg .- vp             # Residual
    vnorm    = sqrt(vg'*(B.*vg))
    vg       = vg./vnorm
    H[i+1,i] = vnorm          
    V[:,i+1] = vg
    v        = reshape(vg,lx1,nel)
  end
  Hb  = H[1:Nev,1:Nev];
  F   = eigen(Hb)

  λ         = F.values
  evec      = F.vectors

# Restart Here


  ϕ         = V[:,1:Nev]*evec

 
else

  λ, ϕ = eigs(Lnew, nev=Nev, ncv=Ncv, which=:LR, maxiter=500)

end  

  Lesshafft_λ = 1.0*im*λ
  
  λi = imag(Lesshafft_λ)
  λr = real(Lesshafft_λ)
  
  p1 = ax1.plot(λr,λi, linestyle="none",marker=".",markersize=16)
  
  #ax1.set_xlim(-6.0, 6.0)
  #ax1.set_ylim(-10.0, 1.0)
  #grid(true)
  
  h2 = figure(num=2,figsize=[8.,6.]);
  ax2 = gca()
  for i in 1:1
   if (manualarnoldi)
     p2 = ax2.plot(xglob,real.(ϕ[:,Nev-i+1]), linestyle="-")
   else   
     p2 = ax2.plot(xglob[r1:r2],real.(ϕ[:,i]), linestyle="-")
   end 

  #  p3 = ax2.plot(xglob[r1:r2],imag.(F.vectors[:,r2-i]), linestyle="-")
  end  



println("Done")










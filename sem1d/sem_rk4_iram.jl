println("Main interface for Time Stepper")
println("Using Implicitly Restarted Arnoldi")

using PolynomialBases
using PyPlot,PyCall
using LinearAlgebra
using IterativeSolvers
using SpecialFunctions
using Roots

include("ArnUpd.jl")

close("all")

# Include the function files
# include("sem_main.jl")

# Local Matrices constructed in Sem_main.jl
# Global Matrices also constructed in Sem_main.jl

# Cg    = QT*Conv*Q    # Global Convection matrix
# Lg    = QT*Lap*Q     # Global Laplacian matrix
# Sg    = QT*Src*Q     # Global Src matrix
# Fg    = QT*Fd*Q      # Global Feedback matrix
# Bg    = QT*B         # Global Mass vector
# Big   = 1.0./Bg      # Global inverse Mass vector
# 
# Oper  = similar(Bg)

rcParams = PyPlot.PyDict(PyPlot.matplotlib."rcParams")

# Analytical Eigenvalues
ω1 = find_zero(airyai,(-3.0,-0.0))
ω2 = find_zero(airyai,(-5.0,-3.0))
ω3 = find_zero(airyai,(-6.0,-5.0))
ω4 = find_zero(airyai,(-7.0,-6.0))
ω5 = find_zero(airyai,(-8.0,-7.0))
ω6 = find_zero(airyai,(-9.5,-8.0))
ω7 = find_zero(airyai,(-10.5,-9.5))

ω  = [ω1, ω2, ω3, ω4, ω5, ω6, ω7]
U  = 6.0
γ  = 1.0 - im*1.0

Ω  = im*(U*U/8 .- U*U/4/γ .+ γ^(1.0/3.0)*(U^(4.0/3.0))/(160.0^(2.0/3.0))*ω)

rcParams["markers.fillstyle"] = "none"
hλ = figure(num=1,figsize=[8.,6.]);
ax1 = gca()
pΛ = plot(real.(Ω),imag.(Ω),linestyle="none",marker="o",markersize=8)

xg    = QT*(vimult.*Geom.xm1[:])

Nev   = 2               # Number of eigenvalues to calculate
EKryl = 2*Nev           # Additional size of Krylov space
LKryl = Nev + EKryl     # Total Size of Krylov space    

V     = zeros(Complex,ndof,LKryl)
H     = zeros(Complex,LKryl,LKryl)

v     = rand(Float64,ndof) + im*rand(Float64,ndof);

# Orthogonalize
# α           = sqrt(v'*(Bg.*v))
# v           = v/α
# V[:,1]      = v
# nkryl       = 1

verbose = true
verbosestep = 100
reortho = 50
ngs     = 2       # Number of Gram-Schmidt
nkryl   = 0
θ       = 0.      # Arnoldi Residual norm

V,H,v,θ,nkryl = ArnUpd!(V,H,Bg,θ,nkryl,LKryl,v,ngs)

cm    = get_cmap("tab10");
rgba0 = cm(0) 
rgba1 = cm(1) 
rgba2 = cm(2) 

dt = 0.0001
plotupd = 200
eigcal  = 500

λn = zeros(Complex,nkryl)

bc = zeros(Complex,1,ndof);
Rhs = similar(v)

rcParams["markers.fillstyle"] = "full"

ifconv = false
Istep  = 0
t = 0.            # Time
i = 0             # Istep

while ~ifconv
  global V,H,v,θ
  global t, i
  global plr,pli
  global OPg
  global nkryl
  global QQ, QRj
  global hλ, ax1, λ, pλ, Ar
  global ifconv

  i = i + 1

  t = t + dt;

#  v = V[:,ik]  

  if verbose && mod(i,verbosestep)==0
    println("Istep=$i, Time=$t")
  end  

# Build the operator only the first time
  if (i==1)
    OPg       = Cg .+ Sg .+ Lg
    for j in 1:ndof
      OPg[j,:] = OPg[j,:]./Bg[j]
    end  
    OPg[1,:] = bc
    OPg[1,1] = 1.0 + im*0.0        # Change operator for BC
  end  

# Apply BC       
  v[1]      = 0.0 + im*0.0
  
# RK4 steps
  v1 = v .+ dt/2.0*OPg*v
  v2 = v .+ dt/2.0*OPg*v1
  v3 = v .+ dt*OPg*v2
  v4 = v .+ dt/6.0*(OPg*(v .+ 2.0*v1 .+ 2.0*v2 .+ v3))

  v  = v4

# Expand Krylov space
  if mod(i,reortho)==0
#   Update Arnoldi Vector 
    V,H,v,θ,nkryl = ArnUpd!(V,H,Bg,θ,nkryl,LKryl,v,ngs)      



#   Perform implicit restart      
    if nkryl == LKryl+1
      global Hb
      kk = nkryl-1
#      Hb = H[1:kk,1:kk]
      F  = eigen(H)          # Uses Lapack routine (dgeev/zgeev)
      fr = real.(F.values)
      fr_sort_i = sortperm(fr,rev=false)   # Increasing order
      μ         = F.values[fr_sort_perm[1:EKryl]]
      nμ        = length(μ)

#      II = Matrix{Complex}(1.0I,kk,kk)
#      QQ = deepcopy(II)
#      for j in 1:1 #EKryl
#        k = fr_sort_i[j]
#        QRj = qr(Hb - F.values[k]*II)          # Expensive. Should replace by Bulge Chase
#        Hb  = QRj.T*QRj.factors + F.values[k]*II
#        QQ  = QQ*QRj.factors
#      end

      QQ,Hs = ShiftedQR(H,μ,nμ,ngs)
      v2    = V*Q[:,Nev+1]        # New Vector
      V     = V*QQ[:,1:Nev]       # Updated Krylov space
      β     = Hs[Nev+1,Nev]       # e_k+1^T*H*e_k
      σ     = QQ[kk,Nev]          # e_k+p^T*Q*e_k
      r     = β*v2 + σ*v          # new residual vector

      H     = Hs

      θ     = sqrt(r'*(Bg.*r))
      v     = r/θ

      nkryl = Nev

      ifconv = true
    end     # nkryl == LKryl+1 

#  end       # mod(i,reortho)    
    if (nkryl == Nev && i>LKryl*reortho)
      if i==reortho*LKryl
#        hλ = figure(num=1,figsize=[8.,6.]);
#        ax1 = gca()
      else
#        ax1.clear()
        pλ[1].remove()
      end  
   
      Ar = V[:,1:Nev]'*(Cg .+ Sg .+ Fg .+ Lg)*V[:,1:Nev]       # V'AV
      λ = eigvals(Ar)
      
      Lesshafft_λ = 1.0*im*λ
      for j in length(λ):-1:1  
        display("$(Lesshafft_λ[j])")
      end  

      pλ = ax1.plot(real.(Lesshafft_λ),imag.(Lesshafft_λ), linestyle="none",marker=".", markersize=8)

      fac = 0.2
      o1 = minimum(real.(Ω))
      l1 = minimum(real.(Lesshafft_λ))
      x1 = min(l1,o1)

      o2 = maximum(real.(Ω))
      l2 = maximum(real.(Lesshafft_λ))
      x2 = max(l2,o2)

      dx = (x2-x1)*fac
      x1 = x1 - dx
      x2 = x2 + dx

      o1 = minimum(imag.(Ω))
      l1 = minimum(imag.(Lesshafft_λ))
      y1 = min(l1,o1) 

      o2 = maximum(imag.(Ω))
      l2 = maximum(imag.(Lesshafft_λ))
      y2 = max(l2,o2)

      dy = (y2-y1)*fac
      y1 = y1 - dy
      y2 = y2 + dy


      ax1.set_xlim((x1,x2))  
      ax1.set_ylim((y1,y2))  

   
      pause(0.001)
    end     # nkryl == Nev  
  end       # mod(i,reortho)    

  end       # mod(i,reortho)    

# Plot
  if mod(i,plotupd)==0
    global hev, ax2
    if i==plotupd
      hev = figure(num=2,figsize=[8.,6.]);
      ax2 = gca()
    end  

    if (i>plotupd)
#       plr[1].remove();
#       pli[1].remove();
#      lo = ax2.get_lines()
      for lo in ax2.get_lines()
        lo.remove()
      end  
    end   

    for ik in 1:nkryl
      plr = ax2.plot(xg,real.(V[:,ik]),color=cm(ik-1));
      pli = ax2.plot(xg,imag.(V[:,ik]),color=cm(ik-1),linestyle="--");
    end  

    pause(0.0001)
  end
 















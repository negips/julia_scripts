println("Main interface for Time Stepper")
println("Using Implicitly Restarted Arnoldi")

using PolynomialBases
using PyPlot,PyCall
using LinearAlgebra
using IterativeSolvers
using SpecialFunctions
using Roots
using Random


include("ArnUpd.jl")
include("ArnIRst.jl")
include("ExplicitShiftedQR.jl")
include("IRAM.jl")

close("all")

# Include the function files
include("sem_main.jl")

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

rng = MersenneTwister(1234)


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

Nev   = 4               # Number of eigenvalues to calculate
EKryl = 3*Nev           # Additional size of Krylov space
LKryl = Nev + EKryl     # Total Size of Krylov space    

vt    = ComplexF64
#vt    = Float64

V     = zeros(vt,ndof,LKryl+1)
Vold  = zeros(vt,ndof,LKryl+1)

H     = zeros(vt,LKryl+1,LKryl)
Hold  = zeros(vt,LKryl+1,LKryl)

v     = randn(vt,ndof);

# Orthogonalize
# α           = sqrt(v'*(Bg.*v))
# v           = v/α
# V[:,1]      = v
# nkryl       = 1

ifarnoldi = true
verbose = false
reortho = 200
verbosestep = reortho #500
nsteps      = 100000

ngs     = 2       # Number of Gram-Schmidt
nkryl   = 0
tol     = 1.0e-12

h,θ,r  = ArnUpd(V,Bg,v,nkryl,ngs)
V[:,1] = r
nkryl  = 1

cm    = get_cmap("tab10");
rgba0 = cm(0) 
rgba1 = cm(1) 
rgba2 = cm(2) 

dt = 0.0001
plotupd = 200

λn = zeros(vt,nkryl)

bc = zeros(vt,1,ndof);
Rhs = similar(v)

rcParams["markers.fillstyle"] = "full"

ifconv = false
Istep  = 0
t = 0.            # Time
i = 0             # Istep

maxouter_it = 50
major_it    = 1

while (~ifconv)
  global V,H,v
  global t, i
  global plr,pli
  global OPg
  global nkryl
  global hλ, ax1, λ, pλ, Ar
  global ifconv
  global major_it
  global Vold,Hold

  local h,r,β
  local U

  i = i + 1

  t = t + dt;

#  v = V[:,ik]  

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

  if (ifarnoldi)
#   Expand Krylov space
    if mod(i,reortho)==0
 #     Update Arnoldi Vector 
 #      V,H,v,θ,nkryl = ArnUpd!(V,H,Bg,θ,nkryl,LKryl,v,ngs)      
       h,β,r = ArnUpd(V,Bg,v,nkryl,ngs)
       H[1:nkryl,nkryl] = h
       H[nkryl+1,nkryl] = β
       nkryl = nkryl + 1
       V[:,nkryl] = r
       v          = r
 
 #     Perform implicit restart

#       V,H,nkryl,β,major_it = IRAM!(V,H,Bg,v,nkryl,LKryl,major_it,Nev,ngs)

       if nkryl == LKryl+1
         Hold = H
         Vold = V
         U,G,nkryl,ifconv = ArnIRst(V,H,Bg,nkryl,LKryl+1,Nev,ngs)
         V = U
         H = G
         
         v = V[:,nkryl]
         β = abs(H[Nev+1,Nev])
         println("Outer Iteration: $major_it; β=$β")
 
         if (β < tol)
           break
         end  

        major_it = major_it + 1
        if (major_it>maxouter_it)
          break
        end  
      end     # nkryl == LKryl+1
      

      if (verbose)
        println("Krylov Size: $nkryl")
      end  

    end       # mod(i,reortho)
  
  elseif i==nsteps
   if verbose && mod(i,verbosestep)==0
    println("Istep=$i, Time=$t")
  end  

   
    break
  end         # ifarnoldi 

end       # while ... 

if (ifarnoldi)
  Hr = H[1:Nev,1:Nev]
  F  = eigen(Hr)
  DT = dt*reortho 
  
  λr = log.(abs.(F.values))/DT
  λi = atan.(imag(F.values),real.(F.values))/DT
  
  λ  = λr .+ im*λi
  
  Lesshafft_λ = 1.0*im*λ
  
  pλ = ax1.plot(real.(Lesshafft_λ),imag.(Lesshafft_λ), linestyle="none",marker=".", markersize=8)
  
  eigvec = V[:,1:Nev]*F.vectors
  
  hv = figure(num=2,figsize=[8.,6.]);
  ax2 = gca()
  for j in 1:Nev
    pvec = ax2.plot(xg,real.(eigvec[:,j]),linestyle="-")
  end  
else
  hv = figure(num=2,figsize=[8.,6.]);
  ax2 = gca()
  pvec = ax2.plot(xg,real.(v),linestyle="-")
end
















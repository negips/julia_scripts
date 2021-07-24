println("Main interface for Time Stepper")
println("Using Implicitly Restarted Arnoldi")

using PyPlot,PyCall
using LinearAlgebra
using IterativeSolvers
using SpecialFunctions
using Roots
using Random
using JLD2


include("ArnUpd.jl")
include("ArnIRst.jl")
include("ExplicitShiftedQR.jl")
include("IRAM.jl")
include("RK4.jl")

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

rng = MersenneTwister(1235)

rcParams = PyPlot.PyDict(PyPlot.matplotlib."rcParams")

# Analytical Eigenvalues
ω1 = find_zero(airyai,(-3.0,-0.0))
ω2 = find_zero(airyai,(-5.0,-3.0))
ω3 = find_zero(airyai,(-6.0,-5.0))
ω4 = find_zero(airyai,(-7.0,-6.0))
ω5 = find_zero(airyai,(-8.0,-7.0))
ω6 = find_zero(airyai,(-9.5,-8.0))
ω7 = find_zero(airyai,(-10.5,-9.5))
ω8 = find_zero(airyai,(-11.8,-10.5))
ω9 = find_zero(airyai,(-12.0,-11.8))
ω10 = find_zero(airyai,(-12.9,-12.0))
ω11 = find_zero(airyai,(-13.8,-12.9))
ω12 = find_zero(airyai,(-14.8,-13.8))
ω13 = find_zero(airyai,(-15.8,-14.8))
ω14 = find_zero(airyai,(-16.8,-15.8))
ω15 = find_zero(airyai,(-17.5,-16.8))

ω  = [ω1, ω2, ω3, ω4, ω5, ω6, ω7, ω8, ω9, ω10, ω11, ω12, ω13, ω14, ω15]
U  = 6.0
γ  = 1.0 - im*1.0

Ω  = im*(U*U/8.0 .- U*U/(4.0*γ) .+ γ^(1.0/3.0)*(U^(4.0/3.0))/(160.0^(2.0/3.0))*ω)

rcParams["markers.fillstyle"] = "none"
hλ = figure(num=1,figsize=[8.,6.]);
ax1 = gca()
pΛ = plot(real.(Ω),imag.(Ω),linestyle="none",marker="o",markersize=8)

xg    = QT*(vimult.*Geom.xm1[:])

Nev   = 20               # Number of eigenvalues to calculate
EKryl = Int64(floor(2.0*Nev))           # Additional size of Krylov space
LKryl = Nev + EKryl     # Total Size of Krylov space    

vt    = ComplexF64
#vt    = Float64

V     = zeros(vt,ndof,LKryl+1)
Vold  = zeros(vt,ndof,LKryl+1)

H     = zeros(vt,LKryl+1,LKryl)
Hold  = zeros(vt,LKryl+1,LKryl)

r     = randn(vt,ndof);

ifarnoldi   = true
ifplot      = false
verbose     = true
reortho     = 2000
verbosestep = reortho #500
nsteps      = 100000
ifsave      = true

ngs     = 0       # Number of Gram-Schmidt
nkryl   = 0
tol     = 1.0e-11

h,θ,v  = ArnUpd(V,Bg,r,nkryl,ngs)
V[:,1] = v
nkryl  = 0

cm    = get_cmap("tab10");
rgba0 = cm(0) 
rgba1 = cm(1) 
rgba2 = cm(2) 

dt = 0.00005

λn = zeros(vt,nkryl)

bc = zeros(vt,1,ndof);
Rhs = similar(v)

rcParams["markers.fillstyle"] = "full"

ifconv = false
t = 0.            # Time
i = 0             # Istep

maxouter_it = 5000
major_it    = 1

if (ifplot)
  hv  = figure(num=2,figsize=[8.,6.]);
  ax2 = gca()
end

# Start iterations
println("Starting Iterations")

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
  global ax2 

  local h,r,β
  local U

  β = 1.0
  i = i + 1

  t = t + dt;

# Build the operator only the first time
  if (i==1)
    OPg       = Cg .+ Sg .+ Lg .+ Fg
    for j in 1:ndof
      OPg[j,:] = OPg[j,:]./Bg[j]
    end  
    OPg[1,:] = bc
    OPg[1,1] = 1.0 + im*0.0        # Change operator for BC
  end  

# Apply BC       
  v[1]      = 0.0 + im*0.0
  
## RK4 steps
#  v1 = v .+ dt/2.0*OPg*v
#  v2 = v .+ dt/2.0*OPg*v1
#  v3 = v .+ dt*OPg*v2
#  v4 = v .+ dt/6.0*(OPg*(v .+ 2.0*v1 .+ 2.0*v2 .+ v3))
#
#  v  = v4
   v  = RK4!(OPg,v,dt)

  if (ifarnoldi)
#   Expand Krylov space
    if mod(i,reortho)==0

       if ifplot
         if (i>reortho) 
           for lo in ax2.get_lines()
             lo.remove()
           end  
         end  

         pv1 = ax2.plot(xg,real.(v),linestyle="-")
       end  

       V,H,nkryl,β,major_it = IRAM!(V,H,Bg,v,nkryl,LKryl,major_it,Nev,ngs)

       v   = V[:,nkryl]

       if (ifplot)
         pv2 = ax2.plot(xg,real.(v),linestyle="--")

         vmin = 1.5*minimum(real.(v))
         vmax = 1.5*maximum(real.(v))
         ax2.set_ylim((vmin,vmax))
         pause(0.001)
         draw() 
       end  

       if (major_it>maxouter_it)
         break
       end  
       if (verbose)
         println("Krylov Size: $nkryl")
       end
       if (β < tol)
         break
       end  

    end       # mod(i,reortho)
  
  else
    if verbose && mod(i,verbosestep)==0
      println("Istep=$i, Time=$t")
    end
   
    if i==nsteps
      break
    end  
  end    # Arnoldi   

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
  
  hev = figure(num=3,figsize=[8.,6.]);
  ax3 = gca()
  for j in 1:Nev
    local pvec1 = ax3.plot(xg,real.(eigvec[:,j]),linestyle="-")
#    local pvec2 = ax3.plot(xg,imag.(eigvec[:,j]),linestyle="--")
  end  
else
  hev = figure(num=3,figsize=[8.,6.]);
  ax3 = gca()
  pvec = ax3.plot(xg,real.(v),linestyle="-")
end



if (ifsave )
  save("nev20_xe40_c0_tol-6.jld2"; AT,N,Nd,xs,xe,nel,U,γ,Ω,xg,vt,Nev,EKryl,LKryl,reortho,V,H,F,DT,λ,Lesshafft_λ)
end  











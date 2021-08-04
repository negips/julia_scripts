println("Main interface for Time Stepper")
println("Using Implicitly Restarted Arnoldi")

using PyPlot,PyCall
using LinearAlgebra
using IterativeSolvers
using SpecialFunctions
using Roots
using Random
using GenericLinearAlgebra          # For eigvals for BigFloat
using Printf
# using JLD2


include("ArnUpd.jl")
include("ArnIRst.jl")
include("ExplicitShiftedQR.jl")
include("BulgeChase.jl")
include("IRAM.jl")
include("RK4.jl")
include("SEM_RK4.jl")

close("all")

ifglobal = false

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
# Parameters defined earlier
#U  = 6.0
#γ  = 1.0 - im*1.0

Ω  = im*(U*U/8.0 .- U*U/(4.0*γ) .+ γ^(1.0/3.0)*(U^(4.0/3.0))/(160.0^(2.0/3.0))*ω)

rcParams["markers.fillstyle"] = "none"
hλ = figure(num=1,figsize=[8.,6.]);
ax1 = gca()
pΛ = plot(real.(Ω),imag.(Ω),linestyle="none",marker="o",markersize=8)

xg          = QT*(vimult.*Geom.xm1[:])
xall        = Geom.xm1[:]
B           = Geom.bm1[:]
Bdssum      = Q*(QT*B)
Binv        = one./Bdssum

ntot        = npts

Nev         = 15                    # Number of eigenvalues to calculate
EKryl       = Int64(floor(2.5*Nev)) # Additional size of Krylov space
LKryl       = Nev + EKryl           # Total Size of Krylov space    
ngs         = 2                   # Number of Gram-Schmidt
tol         = 1.0e-12

vt    = Complex{prec}
#vt    = Float64

V     = zeros(vt,ntot,LKryl+1)
Vold  = zeros(vt,ntot,LKryl+1)

H     = zeros(vt,LKryl+1,LKryl)
Hold  = zeros(vt,LKryl+1,LKryl)

if prec == BigFloat
  r0   = rand(prec,ntot) + im*rand(prec,ntot);
else
  r0   = randn(vt,ntot);
end  

r0     = (one+one*im)sin.(2*pi*xall)
r0[1]  = zro
r      = vimult.*(Q*QT*r0)

ifarnoldi   = false
ifplot      = true
verbose     = true
eigupd      = true
reortho     = 500
verbosestep = reortho #500
nsteps      = 100000
ifsave      = true

nkryl   = 0
h,θ,v  = ArnUpd(V,B,r,nkryl,ngs)
V[:,1] = v
nkryl  = 1

cm    = get_cmap("tab10");
rgba0 = cm(0) 
rgba1 = cm(1) 
rgba2 = cm(2) 

if prec == BigFloat
  dt = BigFloat(0.0001)
else
  dt = 0.0001
end  

λn = zeros(vt,nkryl)

bc = zeros(vt,1,lx1);
Rhs = similar(v)


ifconv = false
t = 0.0*dt        # Time
i = 0             # Istep

maxouter_it = 100
major_it    = 1

rcParams["markers.fillstyle"] = "full"
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
  global hλ, ax1
  global ifconv
  global major_it
  global Vold,Hold
  global ax2 

  local β
  local pλ
  local Hr,evs,λ,λr,λi,Lesshafft_λ,DT,l0

  β = one
  i = i + 1

  t = t + dt;

# Build the operator only the first time
  if (i==1)
    OP[1,:,1] = bc
    OP[1,1,1] = one + im*zro        # Change operator for BC
    B[1,1]    = one + im*zro
    Binv[1,1] = one + im*zro
  end  

# Apply BC
  v[1]      = zro + im*zro
  
  v  = SEM_RK4!(v,dt,nel,lx1,OP,B,Binv,Q,QT,prec)

  if (ifarnoldi)
#   Expand Krylov space
    if mod(i,reortho)==0

       if ifplot
         if (i>reortho) 
           for lo in ax2.get_lines()
             lo.remove()
           end  
         end  

         pv1 = ax2.plot(xall,real.(v),linestyle="-")
       end  

       if nkryl == LKryl
         Hold = H
         Vold = V
       end
       V,H,nkryl,β,major_it = IRAM!(V,H,B,v,nkryl,LKryl,major_it,Nev,ngs)

       v   = V[:,nkryl]

       if (ifplot)
         pv2 = ax2.plot(xall,real.(v),linestyle="--")

         vmin = 1.5*minimum(real.(v))
         vmax = 1.5*maximum(real.(v))
#         ax2.set_ylim((vmin,vmax))
         ax2.set_ylim((-1.0,1.0))
        
         pause(0.001)
         draw() 
       end  

       if (major_it>maxouter_it)
         break
       end  
       if (verbose)
#         println("Major Iteration: $major_it, Krylov Size: $nkryl, β: $β")
        @printf "Major Iteration: %i, Krylov Size: %i, β: %e\n" major_it nkryl β
       end
       if (β < tol)
         println("β = $β")
         break
       end

       if (eigupd) && nkryl == Nev+1

         l0 = ax1.get_lines()
         for il = 2:length(l0)
           l0[il].remove()
         end  

         Hr = H[1:Nev,1:Nev]

         evs = eigvals(Hr)

         DT = dt*reortho 
         
         λr = log.(abs.(evs))/DT
         λi = atan.(imag(evs),real.(evs))/DT
         
         λ  = λr .+ im*λi
         
         Lesshafft_λ = one*im*λ
         
         pλ = ax1.plot(real.(Lesshafft_λ),imag.(Lesshafft_λ), linestyle="none",marker=".", markersize=8)
         ax1.set_xlim((-2.0,6.0))
         ax1.set_ylim((-7.5,0.5))
            
         draw()
         
         pause(0.001)
       end        # eigupd

    end       # mod(i,reortho)
  
  else
    if verbose && mod(i,verbosestep)==0
      println("Istep=$i, Time=$t")
    end
    if (ifplot && mod(i,reortho)==0)
      if (i>reortho) 
        for lo in ax2.get_lines()
          lo.remove()
        end  
      end  
     
      pv1 = ax2.plot(xall,real.(v),linestyle="-")

      vmin = 1.5*minimum(real.(v))
      vmax = 1.5*maximum(real.(v))
      vvmax = max(abs(vmin),abs(vmax))    
      ax2.set_ylim((-vvmax,vvmax))
#      ax2.set_ylim((-1.0,1.0))
      draw()
      pause(0.0001)
    end 
  
    if i==nsteps
      break
    end  
  end    # Arnoldi   

end       # while ... 

if (ifarnoldi)
  Hr = H[1:Nev,1:Nev]

  if prec == BigFloat
    evs = eigvals(Hr)
  else
    F  = eigen(Hr)
    evs = F.values
  end

  DT = dt*reortho 
  
  λr = log.(abs.(evs))/DT
  λi = atan.(imag(evs),real.(evs))/DT
  
  λ  = λr .+ im*λi
  
  Lesshafft_λ = one*im*λ
  
  pλ = ax1.plot(real.(Lesshafft_λ),imag.(Lesshafft_λ), linestyle="none",marker=".", markersize=8)
  
  if prec != BigFloat
    eigvec = V[:,1:Nev]*F.vectors
    
    hev = figure(num=3,figsize=[8.,6.]);
    ax3 = gca()
    for j in 1:Nev
      local pvec1 = ax3.plot(xall,real.(eigvec[:,j]),linestyle="-")
#      local pvec2 = ax3.plot(xg,imag.(eigvec[:,j]),linestyle="--")
    end
  end  
else
  hev = figure(num=3,figsize=[8.,6.]);
  ax3 = gca()
  pvec = ax3.plot(xall,real.(v),linestyle="-")
end



# if (ifsave )
#   save("nev20_xe40_c0_tol-6.jld2"; VT,N,Nd,xs,xe,nel,U,γ,Ω,xg,vt,Nev,EKryl,LKryl,reortho,V,H,F,DT,λ,Lesshafft_λ);
# end  











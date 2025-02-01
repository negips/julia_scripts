println("Main interface for Time Stepper")
println("Using Implicitly Restarted Arnoldi")

using PyPlot,PyCall
using LinearAlgebra
using IterativeSolvers
using SpecialFunctions
using Roots
using Random
# using GenericLinearAlgebra          # For eigvals for BigFloat
# using DoubleFloats
using Printf
using JLD2


include("ArnUpd.jl")
include("ArnIRst.jl")
include("ExplicitShiftedQR.jl")
include("BulgeChase.jl")
include("IRAM.jl")
#include("RK4.jl")
include("$JULIACOMMON/RK4.jl")

close("all")

# Ifglobal
ifglobal = true

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
#Ω  = im*(U*U/8.0 .- U*U/(4.0*γ) .+ γ^(1.0/3.0)*(U^(4.0/3.0))/(160.0^(2.0/3.0))*ω)

Ω0    = im*1.0
U     = 1.0
ϕ     = -π/4.0
R     = 1.0
γ     = R*exp(im*ϕ)
μx    = U/8.0 
μ0    = Ω0 + (U^2)/(4.0*γ) - ((γ*μx*μx)^(1.0/3.0))*ω1 


#cd = imag(γ)
#μ0 = U*U/8.0
#dμ = -(U*U/8.0)*(1.0/cx0)
#Ω  = im*(μ0 .- U*U/(4.0*γ) .+ (γ*dμ*dμ)^(1.0/3.0)*ω)
Ω  = (μ0 .- U*U/(4.0*γ) .+ (γ*μx*μx)^(1.0/3.0)*ω)

# Include the function files
include("sem_main.jl")

rng = MersenneTwister(1235)


xg    = QT*(vimult.*Geom.xm1[:])

Nev         = 5                           # Number of eigenvalues to calculate
EKryl       = Int64(floor(4*Nev))       # Additional size of Krylov space
LKryl       = Nev + EKryl                 # Total Size of Krylov space    
ngs         = 2                           # Number of Gram-Schmidt
tol         = prec(1.0e-10)

vt    = VT # Complex{prec}
#vt    = Float64

V     = zeros(vt,ndof,LKryl+1)
Vold  = zeros(vt,ndof,LKryl+1)

H     = zeros(vt,LKryl+1,LKryl)
Hold  = zeros(vt,LKryl+1,LKryl)

r = randn(vt,ndof)

#if prec == BigFloat
#  r   = rand(prec,ndof) + im*rand(prec,ndof);
#else
#  r   = randn(vt,ndof);
#end  

r     = (one+one*im)sin.(5*pi*xg[:])
r[1]  = prec(0)

ifarnoldi   = true
ifoptimal   = false     # Calculate optimal responses
ifadjoint   = false     # Superceded by ifoptimal
ifplot      = false 
verbose     = true
eigupd      = true
reortho     = 500
if (ifoptimal)
  arnstep   = reortho*2
else
  arnstep   = reortho
end  
verbosestep = arnstep #500
nsteps      = 50000000
ifsave      = false

if (ifadjoint)
  Ω = conj.(Ω)
end  

nkryl   = 0
block   = 1
h,θ,v  = ArnUpd(V,block,Bg,r,nkryl,ngs)
V[:,1] = v
nkryl  = 0

cm    = get_cmap("tab10");
rgba0 = cm(0) 
rgba1 = cm(1) 
rgba2 = cm(2) 

dt = prec(0.0001)

λn = zeros(vt,nkryl)

bc = zeros(vt,1,ndof);
Rhs = similar(v)


# Plots
#-------------------------------------------------- 
rcParams["markers.fillstyle"] = "none"
hλ = figure(num=1,figsize=[8.,6.]);
ax1 = gca()
pΛ = ax1.plot(imag.(Ω),real.(Ω),linestyle="none",marker="o",markersize=8)
Ωi_max = maximum(abs.(imag.(Ω)))
Ωr_min = minimum(real.(Ω))
Ωr_max = maximum(real.(Ω))

ax1.set_xlim((-2*Ωi_max,2*Ωi_max))
ax1.set_ylim((1.1*Ωr_min,1.0))
rcParams["markers.fillstyle"] = "full"
if (ifplot)
  hv  = figure(num=2,figsize=[8.,6.]);
  ax2 = gca()
end
#-------------------------------------------------- 


ifconv = false
t = zro*dt        # Time
i = 0             # Istep

maxouter_it = 150
major_it    = 1

# Start iterations
ifdirect = true
println("Starting Iterations")

while (~ifconv)
  global V,H,v
  global t, i
  global plr,pli
  global OPg, AOPg
  global nkryl
  global hλ, ax1
  global ifconv
  global major_it
  global Vold,Hold,vold
  global ax2
  global ifdirect

  local β

  local pλ
  local Hr,evs,λ,λr,λi,Lesshafft_λ,DT,l0

  β = one
  i = i + 1

  if (ifoptimal) 
    if (mod(i,reortho)==1)
#     Switch every (reortho) steps
      ifdirect = ~ifdirect
#      println("IfDirect: $ifdirect")
    end        
  elseif ifadjoint
    ifdirect = false
  else
    ifdirect = true
  end  

  t = t + dt;

# Build the operator only the first time
# 
  if (i==1)

#   Adjoint Operator BCs    
    AOPg[1,:] = bc
    AOPg[1,1] = one + im*zro        # Change operator for BC
#    OPg       = Cg .+ Sg .+ Lg .+ Fg
#    for j in 1:ndof
#      OPg[j,:] = OPg[j,:]./Bg[j]
#    end  
#   Direct Operator BCs 
    OPg[1,:] = bc
    OPg[1,1] = one + im*zro        # Change operator for BC

  end  

# Apply BC
  v[1]      = zro + im*zro

  if ifdirect
    v       = RK4!(OPg,v,dt)
  else
    v       = RK4!(AOPg,v,dt)
  end

  if (ifarnoldi)

#   Plotting 
    if ifplot && mod(i,reortho)==0
      if (i>reortho) 
        for lo in ax2.get_lines()
          lo.remove()
        end  
      end  
      pv1 = ax2.plot(xg,real.(v),linestyle="-")
    end  

#   Expand Krylov space
    if mod(i,arnstep)==0

      if nkryl == LKryl
        Hold = H
        Vold = V
        vold = v
      end
      V,H,nkryl,β,major_it = IRAM!(V,H,Bg,v,nkryl,LKryl,major_it,Nev,ngs)

      v   = V[:,nkryl]

      if (major_it>maxouter_it)
        break
      end  

      if (verbose)
        @printf "Major Iteration: %3i/%3i, Krylov Size: %3i/%3i, β: %12e\n" major_it maxouter_it nkryl LKryl β
      end
      if (β < tol)
        @printf "Stopping Iteration, β: %12e\n" β
        break
      end  

#     Update Eigenvalues
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
        Lesshafft_λ = λ
        
        pλ = ax1.plot(imag.(Lesshafft_λ),real.(Lesshafft_λ), linestyle="none",marker=".", markersize=8)
        if ifoptimal
          ax1.autoscale(enable=true,axis="both")
        else
          # ax1.set_xlim((-5.0,8.0))
          # ax1.set_ylim((-7.5,2.5))
          ax1.set_xlim((-2.0*Ωi_max,2*Ωi_max))
          ax1.set_ylim((1.1*Ωr_min,1.0))

        end  
           
        draw()
        pause(0.0001)
      end        # eigupd

    end       # mod(i,arnstep)

#   Plotting      
    if (ifplot && mod(i,reortho)==0)
      pv2 = ax2.plot(xg,real.(v),linestyle="--")

      vmin = 1.5*minimum(real.(v))
      vmax = 1.5*maximum(real.(v))
#      ax2.set_ylim((vmin,vmax))
      ax2.set_ylim((-2.0,2.0))
     
      pause(0.00001)
      draw() 
    end  

 
  else
#   No Arnoldi iteration      
    if verbose && mod(i,verbosestep)==0
      println("Istep=$i, Time=$t")
    end
    if (ifplot && mod(i,reortho)==0)
      if (i>reortho) 
        for lo in ax2.get_lines()
          lo.remove()
        end  
      end  
     
      pv1 = ax2.plot(xg,real.(v),linestyle="-")

      vmin = 2.0*minimum(real.(v))
      vmax = 2.0*maximum(real.(v))
      dv   = abs(vmax-vmin)
#      ax2.set_ylim((vmin,vmax))
      ax2.set_ylim((-dv,dv))
     
    end 
  
    if i==nsteps
      break
    end  
  end    # Arnoldi   

end       # while ... 


# Plotting after convergence
if (ifarnoldi)
  Hr = H[1:Nev,1:Nev]

  if prec != Float64
#   eigvals works for BigFloat
#   But eigvecs does not. So if we want the eigenvectors
#   Need to come back to Float64/ComplexF64
#   Or write my own routine for eigenvector calculation.

    evs  = eigvals(Hr)

    H64  = ComplexF64.(Hr)
    F    = eigen(H64)
  else
    F  = eigen(Hr)
    evs = F.values
  end

  DT = dt*reortho 
  
  λr = log.(abs.(evs))/DT
  λi = atan.(imag(evs),real.(evs))/DT
  
  λ  = λr .+ im*λi
  
  Lesshafft_λ = one*λ

  if (eigupd)
    l0 = ax1.get_lines()
    for il = 2:length(l0)
      l0[il].remove()
    end
  end  

  pλ = ax1.plot(imag.(Lesshafft_λ),real.(Lesshafft_λ), linestyle="none",marker=".", markersize=8)
 
# Eigenvectors  
  eigvec = V[:,1:Nev]*F.vectors
  
  hev = figure(num=3,figsize=[8.,6.]);
  ax3 = gca()
  c_map = get_cmap("tab10")
  for j in 1:Nev
    local pvec1 = ax3.plot(xg,real.(eigvec[:,j]),linestyle="-" ,color=c_map(j-1),label="ϕ$(j)(x)")
    local pvec2 = ax3.plot(xg,imag.(eigvec[:,j]),linestyle="-.",color=c_map(j-1),label="ϕ$(j)(x)")
#    local pveca = ax3.plot(xg,abs.(eigvec[:,j]),linestyle="-")
  end
  legend()

#  Ar        = eigvec'*diagm(Bg)*OPg*eigvec
#  λ_opt     = one*im*eigvals(Ar);
#  pλ2       = ax1.plot(real.(λ_opt),imag.(λ_opt), linestyle="none",marker=".", markersize=8)
#  ax1.set_xlim((-4.0,8.0))
#  ax1.set_ylim((-7.5,2.5))

else
  hev = figure(num=3,figsize=[8.,6.]);
  ax3 = gca()
  pvec = ax3.plot(xg,real.(v),linestyle="-")
end

vnorm = norm(eigvec'*diagm(Bg)*eigvec - I)
@printf("Vnorm: %12e", vnorm)

if (ifsave)
  save("nev15_corr_e-2_double64.jld2","evs",evs, "evec",eigvec);
end  


println("Done.")








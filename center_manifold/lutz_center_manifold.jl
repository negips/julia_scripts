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
using IterativeSolvers
using Peaks
using Statistics


include("ArnUpd.jl")
include("ArnIRst.jl")
include("ExplicitShiftedQR.jl")
include("BulgeChase.jl")
include("IRAM.jl")
#include("RK4.jl")
include("$JULIACOMMON/RK4.jl")
include("OP_RK4.jl")

close("all")

# Ifglobal
ifglobal = true

rcParams = PyPlot.PyDict(PyPlot.matplotlib."rcParams")

# Analytical Eigenvalues
ω1    = find_zero(airyai,(-3.0,-0.0))
ω2    = find_zero(airyai,(-5.0,-3.0))
ω3    = find_zero(airyai,(-6.0,-5.0))
ω4    = find_zero(airyai,(-7.0,-6.0))
ω5    = find_zero(airyai,(-8.0,-7.0))
ω6    = find_zero(airyai,(-9.5,-8.0))
ω7    = find_zero(airyai,(-10.5,-9.5))
ω8    = find_zero(airyai,(-11.8,-10.5))
ω9    = find_zero(airyai,(-12.0,-11.8))
ω10   = find_zero(airyai,(-12.9,-12.0))
ω11   = find_zero(airyai,(-13.8,-12.9))
ω12   = find_zero(airyai,(-14.8,-13.8))
ω13   = find_zero(airyai,(-15.8,-14.8))
ω14   = find_zero(airyai,(-16.8,-15.8))
ω15   = find_zero(airyai,(-17.5,-16.8))

ω     = [ω1, ω2, ω3, ω4, ω5, ω6, ω7, ω8, ω9, ω10, ω11, ω12, ω13, ω14, ω15]
#Ω  = im*(U*U/8.0 .- U*U/(4.0*γ) .+ γ^(1.0/3.0)*(U^(4.0/3.0))/(160.0^(2.0/3.0))*ω)

Ω0    = im*1.0
U     = 1.0
ϕ     = -π/4.0
R     = 1.0
γ     = R*exp(im*ϕ)
μx    = U/8.0 
μ0    = Ω0 + (U^2)/(4.0*γ) - ((γ*μx*μx)^(1.0/3.0))*ω1 
δ5    = (-1.0 + 1.0im)*0.1

ifconj = false
if (ifconj)
  U    = conj(U)
  γ    = conj(γ)
  μ0   = conj(μ0)
  μx   = conj(μx)
end  

#cd = imag(γ)
#μ0 = U*U/8.0
#dμ = -(U*U/8.0)*(1.0/cx0)
#Ω  = im*(μ0 .- U*U/(4.0*γ) .+ (γ*dμ*dμ)^(1.0/3.0)*ω)
Ω  = (μ0 .- U*U/(4.0*γ) .+ (γ*μx*μx)^(1.0/3.0)*ω)

# Include the function files
include("sem_main.jl")

rng = MersenneTwister(1235)


xg          = QT*(vimult.*Geom.xm1[:])
Nev         = 5                           # Number of eigenvalues to calculate
EKryl       = Int64(floor(4*Nev))         # Additional size of Krylov space
LKryl       = Nev + EKryl                 # Total Size of Krylov space    
ngs         = 2                           # Number of Gram-Schmidt
tol         = prec(1.0e-10)

vt          = VT # Complex{prec}
#vt         = Float64

ndofE       = ndof + 1
V           = zeros(vt,ndofE,LKryl+1)
Vold        = zeros(vt,ndofE,LKryl+1)

H           = zeros(vt,LKryl+1,LKryl)
Hold        = zeros(vt,LKryl+1,LKryl)


#r           = (one+one*im)sin.(5*pi*xg[:])
#r[1]        = prec(0)

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
ifsave      = true

if (ifadjoint)
  Ω = conj.(Ω)
end  


## Plots
##-------------------------------------------------- 
#rcParams["markers.fillstyle"] = "none"
#hλ = figure(num=1,figsize=[8.,6.]);
#ax1 = gca()
#pΛ = ax1.plot(imag.(Ω),real.(Ω),linestyle="none",marker="o",markersize=8)
#Ωi_max = maximum(abs.(imag.(Ω)))
#Ωr_min = minimum(real.(Ω))
#Ωr_max = maximum(real.(Ω))
#
#ax1.set_xlim((-2*Ωi_max,2*Ωi_max))
#ax1.set_ylim((1.1*Ωr_min,1.0))
#rcParams["markers.fillstyle"] = "full"
#if (ifplot)
#  hv  = figure(num=2,figsize=[8.,6.]);
#  ax2 = gca()
#end
##-------------------------------------------------- 



Dir         = load("direct_GL_nev5.jld2")
Dir_conj    = load("direct_conj_GL_nev5.jld2")

Adj         = load("adjoint_GL_nev5.jld2")
Adj_conj    = load("adjoint_conj_GL_nev5.jld2")

Y15         = load("direct_resolvent_GL_y15_nev5.jld2")
Y25         = load("direct_resolvent_GL_y25_nev5.jld2")


λ           = Dir["evs"]
i1          = argmax(real.(λ))
v1          = Dir["evec"][:,i1]
λ1          = λ[i1]

λ           = Dir_conj["evs"]
i2          = argmax(real.(λ))
v2          = Dir_conj["evec"][:,i2]
λ2          = λ[i2]

λa          = Adj["evs"]
ia1         = argmax(real.(λa))
w1          = Adj["evec"][:,ia1]
λa1         = λa[ia1]

λa          = Adj_conj["evs"]
ia2         = argmax(real.(λa))
w2          = Adj_conj["evec"][:,ia2]
λa2         = λa[ia2]

λ           = Y15["evs"]
i1          = argmax(real.(λ))
y15E        = Y15["evec"][:,i1]
y15         = y15E[1:ndof]/y15E[ndofE]

λ           = Y25["evs"]
i1          = argmax(real.(λ))
y25E        = Y25["evec"][:,i1]
y25         = y25E[1:ndof]/y25E[ndofE]

#---------------------------------------------------------------------- 
xg          = QT*(vimult.*Geom.xm1[:])

# Normalize vectors
α1    = v1'*(Bg.*v1)
v1    = v1./α1
α1    = (w1'*(Bg.*v1))'
w1    = w1./α1

# α2    = v2'*(Bg.*v2)
# v2    = v2./α2
# α2    = (w2'*(Bg.*v2))'
# w2    = w2./α2


# Plotting
#---------------------------------------------------------------------- 
close("all")
hev         = figure(num=3,figsize=[8.,6.]);
ax3         = gca()
c_map       = get_cmap("tab10")
pv1         = ax3.plot(xg,real.(v1),linestyle="-",color=c_map(0),label=L"\mathfrak{R}(ϕ_{c})")
#pv2         = ax3.plot(xg,imag.(v1),linestyle="-.",color=c_map(0),label=L"\mathfrak{Im}(v)") #
pv2         = ax3.plot(xg,imag.(v1),linestyle="--",color=c_map(0)) #

#pv3         = ax3.plot(xg,real.(v2),linestyle="--",color=c_map(1),label="V2(x)")
#pv4         = ax3.plot(xg,imag.(v2),linestyle="-.",color=c_map(1)) # ,label="ϕ$(j)(x)")

pw1         = ax3.plot(xg,real.(w1),linestyle="-",color=c_map(1),label=L"\mathfrak{R}(χ_{c})")
#pw2         = ax3.plot(xg,imag.(w1),linestyle="-.",color=c_map(1),label=L"\mathfrak{Im}(w)") #
pw2         = ax3.plot(xg,imag.(w1),linestyle="--",color=c_map(1)) #
 
# pw3         = ax3.plot(xg,real.(w2),linestyle="--",color=c_map(3),label="W2(x)")
# pw4         = ax3.plot(xg,imag.(w2),linestyle="-.",color=c_map(3)) # ,label="ϕ$(j)(x)")

py15_1      = ax3.plot(xg,real.(y15),linestyle="-",color=c_map(2),label=L"\mathfrak{R}(y^{1}_{15})")
#py15_2      = ax3.plot(xg,imag.(y15),linestyle="-.",color=c_map(2),label=L"\mathfrak{Im}(y_{15})") #
py15_2      = ax3.plot(xg,imag.(y15),linestyle="--",color=c_map(2)) #

py25_1      = ax3.plot(xg,real.(y25),linestyle="-",color=c_map(3),label=L"\mathfrak{R}(y^{1}_{25})")
#py25_2      = ax3.plot(xg,imag.(y25),linestyle="-.",color=c_map(3),label=L"\mathfrak{Im}(y_{25})") #
py25_2      = ax3.plot(xg,imag.(y25),linestyle="--",color=c_map(3)) #

ax3.set_xlabel(L"x",fontsize=16)
ax3.set_ylabel(L"A",fontsize=16)

legend(fontsize=14)
hev.savefig("GL_fields.eps")
#---------------------------------------------------------------------- 

vl          = Q*v1
wl          = Q*w1
y15l        = Q*y15
y25l        = Q*y25
xl          = Q*xg

z15         = (wl)'*(B.*(xl.*vl))
z115        = (wl)'*(B.*(xl.*y15l))
z125        = (wl)'*(B.*(xl.*y25l))

z25         = (wl)'*(B.*(SLap*(vl)))
z215        = (wl)'*(B.*(SLap*(y15l)))
z225        = (wl)'*(B.*(SLap*(y25l)))

z555        = (wl)'*(B.*((conj.(vl).*vl).*vl))

@printf("z_{15}  : %5.3f %5.3fi\n",real(z15) , imag(z15))
@printf("z_{115} : %5.3f %5.3fi\n",real(z115), imag(z115))
@printf("z_{125} : %5.3f %5.3fi\n",real(z125), imag(z125))

@printf("z_{25}  : %5.3f %5.3fi\n",real(z25) , imag(z25))
@printf("z_{215} : %5.3f %5.3fi\n",real(z215), imag(z215))
@printf("z_{225} : %5.3f %5.3fi\n",real(z225), imag(z225))

@printf("z_{555} : %5.3f %5.3fi\n",real(z555), imag(z555))

# Spectra Plot

# Plots
#-------------------------------------------------- 
rcParams["markers.fillstyle"] = "none"
hλ          = figure(num=1,figsize=[6.0,5.]);
ax1         = gca()
nω          = 10
pΛ1         = ax1.plot(imag.(Ω[1:nω]),real.(Ω[1:nω]),linestyle="none",marker="o",markersize=8,label="Analytical - [0,∞)")
Ωi_max      = maximum(abs.(imag.(Ω[1:nω])))
Ωr_min      = minimum(real.(Ω[1:nω]))
Ωr_max      = maximum(real.(Ω[1:nω]))

ax1.set_xlim((0.0,1.1*Ωi_max))
ax1.set_ylim((1.1*Ωr_min,0.5))
rcParams["markers.fillstyle"] = "full"

# Arnoldi Eigenvalues
λ           = Dir["evs"]
pΛ2         = ax1.plot(imag.(λ),real.(λ),linestyle="none",marker="o",markersize=4,label="Numerical - [0,40]")
ax1.set_xlabel(L"\mathfrak{Im}(ω)",fontsize=18)
ax1.set_ylabel(L"\mathfrak{R}(ω)",fontsize=18)
pleg = legend(loc="center left",fontsize=12)
hλ.savefig("GL_spectra.eps")


#---------------------------------------------------------------------- 
function StuartLandau(Z::T,ω,μ1,μ2,δ5,z1,z11,z12,z2,z22,z21,zzz) where {T <: Number}

  DZ = T(0)
  DZ = ω*Z + μ1*z1*Z + (μ1^2)*z11*Z + μ1*μ2*(z12 + z21)*Z + μ2*z2*Z + (μ2^2)*z22*Z + δ5*zzz*conj.(Z).*Z.*Z

  return DZ
end
#---------------------------------------------------------------------- 

ωc           = Ω[1]
SLμ(x,μ1,μ2) = StuartLandau(x,ωc,μ1,μ2,δ5,z15,z115,z125,z25,z225,z215,z555) 


# Saved values (Ginzburg-Landau)
#-------------------------------------------------- 
δμx_v     = zeros(Float64,8)
Ω_v1      = zeros(Float64,8)

δμx_v[1]  = -0.0000
Ω_v1[1]   =  1.0000

δμx_v[2]  = -0.00313
Ω_v1[2]   =  1.0139

δμx_v[3]  = -0.00625
Ω_v1[3]   =  1.0272

δμx_v[4]  = -0.0125
Ω_v1[4]   =  1.0519

δμx_v[5]  = -0.01875
Ω_v1[5]   =  1.0753

δμx_v[6]  = -0.025
Ω_v1[6]   =  1.0977

δμx_v[7]  = -0.03125
Ω_v1[7]   =  1.1194

δμx_v[8]  = -0.0375
Ω_v1[8]   =  1.1403

δγ_v      = zeros(Float64,11)
Ω_v2      = zeros(Float64,11)

δγ_v[1]   = -0.0000
Ω_v2[1]   =  1.0000

δγ_v[2]   = -0.0250
Ω_v2[2]   =  1.0037

δγ_v[3]   = -0.0500
Ω_v2[3]   =  1.0081

δγ_v[4]   = -0.1000
Ω_v2[4]   =  1.0174

δγ_v[5]   = -0.1500
Ω_v2[5]   =  1.0278

δγ_v[6]   = -0.2000
Ω_v2[6]   =  1.0401

δγ_v[7]   = -0.2500
Ω_v2[7]   =  1.0542

δγ_v[8]   = -0.3000
Ω_v2[8]   =  1.0708

δγ_v[9]   = -0.3500
Ω_v2[9]   =  1.0899

δγ_v[10]  = -0.4000
Ω_v2[10]  =  1.1117

δγ_v[11]  = -0.4500
Ω_v2[11]  =  1.1368
#-------------------------------------------------- 

# μ' variation
@printf("\n μ' variation: δ'_{3} \n")

dt          = 0.001
nsteps      = 1000000
histstep    = 10
nhist       = Int(nsteps/histstep)

zhist       = zeros(vt,nhist)
Time        = zeros(Float64,nhist)

t = 0.0
z = vt(0.01)

SL_Ω_v1     = 0.0*Ω_v1
SL_Ω_v1[1]  = 1.0

for jj in 2:length(δμx_v)

  global z,t
  
  δμ1         = -δμx_v[jj]
  δμ2         = 0.0
  SL(x)       = SLμ(x,δμ1,δμ2)
  
  t = 0.0
  z = vt(0.01)
  for i in 1:nsteps
  
    t = t + dt
    z = OP_RK4(SL,z,dt)
  
    if (mod(i,histstep) == 0)
      j = Int(i/histstep)
      zhist[j] = z
      Time[j] = t
    end  
  
  end

  # Frequency
  peak_ind    = argmaxima(real.(zhist))
  peak_times  = Time[peak_ind]
  delta_times = diff(peak_times)
  afreq       = 2.0*π./delta_times
  npeaks      = length(afreq)
  if npeaks > 9
    ωend      = afreq[npeaks-9:npeaks]
  else
    ωend      = afreq
    println("Not enough cycles: $npeaks")
  end

  ωmean       = mean(ωend)
  @printf("δμx : %.5f ; Ω: %.4e; ΩGL: %.4e \n", δμ1, ωmean, Ω_v1[jj])
  SL_Ω_v1[jj] = ωmean
end

h4  = figure(num=4,figsize=[8.0,7.0])
ax4 = gca() 
ax4.plot(-δμx_v/abs(μx),Ω_v1,linestyle="none",marker="o",markersize=8,label="Ginzburg-Landau")
rcParams["markers.fillstyle"] = "none"
ax4.plot(-δμx_v/abs(μx),SL_Ω_v1,linestyle="none",marker="s",markersize=8,label="Center-Manifold")
rcParams["markers.fillstyle"] = "full"

ax4.set_ylabel(L"Ω",fontsize=18)
ax4.set_xlabel(L"δ'_{3}/|δ^{0}_{3}|",fontsize=18)
pleg4 = legend(loc="upper left",fontsize=12)
h4.savefig("center_manifold_omega1.eps")


# Viscosity variation
@printf("\n γ (viscosity) variation: δ'_{4} \n")
nsteps      = 2000000
nhist       = Int(nsteps/histstep)

zhist       = zeros(vt,nhist)
Time        = zeros(Float64,nhist)

SL_Ω_v2     = 0.0*Ω_v2
SL_Ω_v2[1]  = 1.0

for jj in 2:length(δγ_v)

  global z,t
  
  δμ1         = 0.0
  δμ2         = δγ_v[jj]
  SL(x)       = SLμ(x,δμ1,δμ2)
  
  t = 0.0
  z = vt(0.01)
  for i in 1:nsteps
  
    t = t + dt
    z = OP_RK4(SL,z,dt)
  
    if (mod(i,histstep) == 0)
      j = Int(i/histstep)
      zhist[j] = z
      Time[j] = t
    end  
  
  end

  # Frequency
  peak_ind    = argmaxima(real.(zhist))
  peak_times  = Time[peak_ind]
  delta_times = diff(peak_times)
  afreq       = 2.0*π./delta_times
  npeaks      = length(afreq)
  if npeaks > 9
    ωend      = afreq[npeaks-9:npeaks]
  else
    ωend      = afreq
    println("Not enough cycles: $npeaks")
  end

  ωmean       = mean(ωend)
  @printf("δγ : %.5f ; Ω: %.4e; ΩGL: %.4e \n", -δμ2, ωmean, Ω_v2[jj])
  SL_Ω_v2[jj] = ωmean
end

h5  = figure(num=5,figsize=[8.0,7.0])
ax5 = gca() 
ax5.plot(-δγ_v/abs(γ),Ω_v2,linestyle="none",marker="o",markersize=8,label="Ginzburg-Landau")
rcParams["markers.fillstyle"] = "none"
ax5.plot(-δγ_v/abs(γ),SL_Ω_v2,linestyle="none",marker="s",markersize=8,label="Center-Manifold")
rcParams["markers.fillstyle"] = "full"
ax5.set_ylabel(L"Ω",fontsize=18)
ax5.set_xlabel(L"-δ'_{4}/|δ^{0}_{4}|",fontsize=18)
pleg5 = legend(loc="upper left",fontsize=12)
h5.savefig("center_manifold_omega2.eps")

println("Done.")












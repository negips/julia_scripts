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

α2    = v2'*(Bg.*v2)
v2    = v2./α2
α2    = (w2'*(Bg.*v2))'
w2    = w2./α2


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

py15_1      = ax3.plot(xg,real.(y15),linestyle="-",color=c_map(2),label=L"\mathfrak{R}(y_{15})")
#py15_2      = ax3.plot(xg,imag.(y15),linestyle="-.",color=c_map(2),label=L"\mathfrak{Im}(y_{15})") #
py15_2      = ax3.plot(xg,imag.(y15),linestyle="--",color=c_map(2)) #

py25_1      = ax3.plot(xg,real.(y25),linestyle="-",color=c_map(3),label=L"\mathfrak{R}(y_{25})")
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




println("Done.")












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

ifconj = true
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

vt    = VT # Complex{prec}
#vt    = Float64

V     = zeros(vt,ndof,LKryl+1)
Vold  = zeros(vt,ndof,LKryl+1)

H     = zeros(vt,LKryl+1,LKryl)
Hold  = zeros(vt,LKryl+1,LKryl)


r     = (one+one*im)sin.(5*pi*xg[:])
r[1]  = prec(0)

ifarnoldi   = true
ifoptimal   = false     # Calculate optimal responses
ifadjoint   = true      # Superceded by ifoptimal
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


# # Plotting after convergence
# if (ifarnoldi)
#   Hr = H[1:Nev,1:Nev]
# 
#   if prec != Float64
# #   eigvals works for BigFloat
# #   But eigvecs does not. So if we want the eigenvectors
# #   Need to come back to Float64/ComplexF64
# #   Or write my own routine for eigenvector calculation.
# 
#     evs  = eigvals(Hr)
# 
#     H64  = ComplexF64.(Hr)
#     F    = eigen(H64)
#   else
#     F  = eigen(Hr)
#     evs = F.values
#   end
# 
#   DT = dt*reortho 
#   
#   λr = log.(abs.(evs))/DT
#   λi = atan.(imag(evs),real.(evs))/DT
#   
#   λ  = λr .+ im*λi
#   
#   Lesshafft_λ = one*λ
# 
#   if (eigupd)
#     l0 = ax1.get_lines()
#     for il = 2:length(l0)
#       l0[il].remove()
#     end
#   end  
# 
#   pλ = ax1.plot(imag.(Lesshafft_λ),real.(Lesshafft_λ), linestyle="none",marker=".", markersize=8)
#  
# # Eigenvectors  
#   eigvec = V[:,1:Nev]*F.vectors
#   
#   hev = figure(num=3,figsize=[8.,6.]);
#   ax3 = gca()
#   c_map = get_cmap("tab10")
#   for j in 1:Nev
#     local pvec1 = ax3.plot(xg,real.(eigvec[:,j]),linestyle="-" ,color=c_map(j-1),label="ϕ$(j)(x)")
#     local pvec2 = ax3.plot(xg,imag.(eigvec[:,j]),linestyle="-.",color=c_map(j-1)) # ,label="ϕ$(j)(x)")
# #    local pveca = ax3.plot(xg,abs.(eigvec[:,j]),linestyle="-")
#   end
#   legend()
# 
# #  Ar        = eigvec'*diagm(Bg)*OPg*eigvec
# #  λ_opt     = one*im*eigvals(Ar);
# #  pλ2       = ax1.plot(real.(λ_opt),imag.(λ_opt), linestyle="none",marker=".", markersize=8)
# #  ax1.set_xlim((-4.0,8.0))
# #  ax1.set_ylim((-7.5,2.5))
# 
# else
#   hev = figure(num=3,figsize=[8.,6.]);
#   ax3 = gca()
#   pvec = ax3.plot(xg,real.(v),linestyle="-")
# end


Dir         = load("direct_GL_nev5.jld2")
Dir_conj    = load("direct_conj_GL_nev5.jld2")

Adj         = load("adjoint_GL_nev5.jld2")
Adj_conj    = load("adjoint_conj_GL_nev5.jld2")

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

xg          = QT*(vimult.*Geom.xm1[:])

close("all")
hev         = figure(num=3,figsize=[8.,6.]);
ax3         = gca()
c_map       = get_cmap("tab10")
pv1         = ax3.plot(xg,real.(v1),linestyle="-",color=c_map(0),label=L"\mathfrak{R}(v)")
pv2         = ax3.plot(xg,imag.(v1),linestyle="-.",color=c_map(0),label=L"\mathfrak{Im}(v)") # ,label="ϕ$(j)(x)")

#pv3         = ax3.plot(xg,real.(v2),linestyle="--",color=c_map(1),label="V2(x)")
#pv4         = ax3.plot(xg,imag.(v2),linestyle="-.",color=c_map(1)) # ,label="ϕ$(j)(x)")

pw1         = ax3.plot(xg,real.(w1),linestyle="-",color=c_map(1),label=L"\mathfrak{R}(w)")
pw2         = ax3.plot(xg,imag.(w1),linestyle="-.",color=c_map(1),label=L"\mathfrak{Im}(w)") # ,label="ϕ$(j)(x)")
# 
# pw3         = ax3.plot(xg,real.(w2),linestyle="--",color=c_map(3),label="W2(x)")
# pw4         = ax3.plot(xg,imag.(w2),linestyle="-.",color=c_map(3)) # ,label="ϕ$(j)(x)")

legend()


# Normalize vectors
α1    = v1'*(Bg.*v1)
v1    = v1./α1
α1    = (w1'*(Bg.*v1))'
w1    = w1./α1

α2    = v2'*(Bg.*v2)
v2    = v2./α2
α2    = (w2'*(Bg.*v2))'
w2    = w2./α2


Q1          = (I - (Bg.*v1)*w1')
Q2          = (I - v2*w2')

#g15_1       = (Lg*v1)./γ
temp        = SLap*(Q*v1)
temp[1]     = 0.0
g15         = temp - (Q*v1)*(Q*w1)'*(B.*temp)
δ15         = w1'*QT*(B.*g15)
g15_1       = QT*(B.*g15)
# QT*B.*(SLap*(Q*v1) -
#g15_1[1]    = 0.0       # Boundary condition

temp        = Q*(xg.*v1)
temp[1]     = 0.0
g25         = temp - (Q*v1)*(Q*w1)'*(B.*temp)
δ25         = w1'*QT*(B.*g25)
g25_1       = QT*(B.*g25)


bc          = zeros(vt,1,ndof);
OPs         = (0.0 + λ1)*I - OPg 
OPs[1,:]    = bc
OPs[1,1]    = one + im*zro        # Change operator for BC
OPs         = Q1*OPs*Q1
y15         = copy(g15_1)
y15         = gmres!(y15,OPs,g15_1,reltol=1.0e-10,restart=30,verbose=true) 

OPs         = (0.0 + λ1)*I - OPg 
OPs[1,:]    = bc
OPs[1,1]    = one + im*zro        # Change operator for BC
OPs         = Q1*OPs*Q1
y25         = copy(g25_1)
y25         = gmres!(y25,OPs,g25_1,reltol=1.0e-10,orth_meth=DGKS(),restart=30,maxiter=1000,verbose=true) 

# Build Forcings


ax3.plot(xg,real.(y15),linestyle="-" ,color=c_map(2),label=L"\mathfrak{R}(y_{15})")
ax3.plot(xg,imag.(y15),linestyle="-.",color=c_map(2),label=L"\mathfrak{Im}(y_{15})")

ax3.plot(xg,real.(y25),linestyle="-" ,color=c_map(3),label=L"\mathfrak{R}(y_{25})")
ax3.plot(xg,imag.(y25),linestyle="-." ,color=c_map(3),label=L"\mathfrak{Im}(y_{25})")

legend()

println("Done.")












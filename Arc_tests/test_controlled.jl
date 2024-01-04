#!/bin/julia

println("Testing Arc branches")

using LinearAlgebra
using Random
using PyPlot

include("$JULIACOMMON/MoveFigure.jl")

close("all")

lafs = 16

#rng = MersenneTwister(1234);
Random.seed!(1238);

N1          = 4
N2          = 30 
N3          = 16         # Non-normal subspace


@assert N1 % 2 == 0
@assert N2 % 2 == 0
@assert N3 % 2 == 0

N1by2       = Int(N1/2)
N2by2       = Int(N2/2)

N           = N1 + N2
Nby2        = Int(N/2)
σ           = [-0.1*rand(N1by2); -rand(N2by2)]
#σ[1]        = -0.001
ω           = 0.1*rand(Nby2)
Ω           = zeros(ComplexF64,N)

A           = zeros(Float64,N,N)
for i in 1:Nby2
  j = (i-1)*2 + 1
  A[j,j]          = σ[i]
  A[j,j+1]        = ω[i]
  A[j+1,j+1]      = σ[i]
  A[j+1,j]        = -ω[i]

  Ω[j]            = σ[i] + im*ω[i]
  Ω[j+1]          = σ[i] - im*ω[i]
end

V           = zeros(Float64,N,N)
W           = zeros(Float64,N,N)

v           = rand(Float64,N)
vnorm       = norm(v)
v           = v/vnorm
V[:,1]      = v
W[:,1]      = v

for i in 2:N
  w         = rand(Float64,N)
  α         = W[:,1:i-1]'*w
  w         = w .- W[:,1:i-1]*α
  wnorm     = norm(w)
  w         = w./wnorm
  W[:,i]    = w
 
  if i<=N3
    θ       = rand(Float64,i)
    local v = W[:,1:i]*θ
    V[:,i]  = v
  else
    V[:,i]  = w
  end  

end  

L            = V*A*inv(V)
He           = norm(L*L' .- L'*L)/norm(L*L)
println("Henrici Index: $He")



t            = 2.0

expL         = exp(t*L)
λ            = eigvals(expL)
λr           = real.(λ)
λi           = imag.(λ)

Ωr           = real.(Ω)
Ωi           = imag.(Ω)

EΩ           = exp.(t*Ω)
EΩr          = real.(EΩ)
EΩi          = imag.(EΩ)

h1           = figure(num=1)
ax1          = h1.subplots()
plλ          = ax1.plot(λi,λr,linestyle="none",marker="o",markersize=6,fillstyle="none")
plΩ          = ax1.plot(EΩi,EΩr,linestyle="none",marker="x",markersize=6,fillstyle="none")

ϕ            = LinRange(0,2.0*π,10000)
cs           = cos.(ϕ)
sn           = sin.(ϕ)
circpl       = ax1.plot(cs,sn,linestyle="--")

MoveFigure(h1,1250,500)


# Optimal tests
#--------------------------------------------------  

nT           = 20
T0           = 1.0
Tend         = 10.0
T            = LinRange(T0,Tend,nT)

h2           = figure(num=2)
ax2          = h2.subplots()

m            = 8              # Reduced Order Projection subspace size
Σ            = zeros(nT,m)
Λ            = zeros(ComplexF64,nT,m)

cm           = get_cmap("tab10") 

for i in 1:nT

  global Σ

  local t      = T[i]
  local expL   = exp(t*L)

  F            = eigen(expL*expL')
  S            = abs.(F.values)
  ind          = sortperm(S,rev=true)
  Σ[i,:]       = S[ind[1:m]]

  vecs         = F.vectors[:,ind[1:m]]
 
  global Lred         = vecs'*L*vecs
  λ2           = eigvals(Lred)
  λa           = abs.(λ2)
  λ2r          = 1.0/t*log.(λa)
  λ3           = λ2./λa
  #λ2i          = (1.0/t)*atan.(imag.(λ3),real.(λ3)) 
  λ2i          = (1.0/t)*atan.(imag.(λ3)./real.(λ3)) 
  Λ[i,:]       = (λ2r .+ im*λ2i)

#  Λ[i,:]       = 1.0/t*log.(Complex.(λ2))
  
#  ax2.plot(sqrt.(abs.(S)),linestyle="none",marker="o",markersize=6,fillstyle="none")
end

for i in 1:m
  ax2.plot(T,sqrt.(Σ[:,i]),linestyle="-",color=cm(i))
end
ax2.set_ylabel(L"\Sigma",fontsize=lafs)
ax2.set_xlabel("T",fontsize=lafs)

#ax2.set_yscale("log")

h3           = figure(num=3)
ax3          = h3.subplots()
pl0          = ax3.plot(Ωi,Ωr,linestyle="none",marker="o",markersize=5,fillstyle="full",color="red")

for i in 1:m
  plΛ        = ax3.plot(imag.(Λ[:,i]),real.(Λ[:,i]),linestyle="none",marker="x",markersize=4,fillstyle="none",color=cm(i))
end
MoveFigure(h3,600,500)


# Compute Pseudo-spectra
nz           = 800
zi           = LinRange(0.0,0.25,nz)
zr           = LinRange(-1.2,0.1,nz)

#zi           = LinRange(0.09,0.11,nz)
#zr           = LinRange(-0.64,-0.56,nz)

pseudospectra  = Matrix{Float64}(undef,nz,nz)

for j in 1:nz
  for i in 1:nz
    zz = zr[i] + im*zi[j]
    F = svd(zz*I - L)
    pseudospectra[i,j] = minimum(F.S)
  end
end  

cont = ax3.contour(zi,zr,log.(pseudospectra),50)
colorbar(cont)

# Numerical Range




println("Done.")









#!/bin/julia
using LinearAlgebra
using Roots
using PyPlot
#using Symbolics

include("NullClines.jl")
include("NullClineFcn.jl")
include("ElementOf.jl")
include("NullClineParams.jl")
include("MoveFigure.jl")

function evs(k,α,ν1,ν2)
  @assert size(α) == (2,2)

  LMat =  [(α[1,1] - ν1*k^2)   α[1,2];
           α[2,1]             (α[2,2] - ν2*k^2)]
  ee   = eigvals(LMat)

  er   = real.(ee)
  ei   = imag.(ee)

  return er,ei
end  

#-------------------------------------------------- 

include("build_nullcline_fcns.jl")
lafs              = 16

# Wavenumber Instability
ϵ = 1.0
η = 1.0

α      = zeros(Float64,2,2)
α[1,1] = pars.gcx[1]*η
α[1,2] = pars.gcy[1]*η

α[2,1] = pars.fcx[1]/ϵ
α[2,2] = pars.fcy[1]/ϵ

ν1  = 1.0
ν2  = 1.0

ωr  = 1.0
ωi  = 1.0
ω   = 1.0

LMat(k,ω) =  [(α[1,1] - ν1*k^2 - ω)          α[1,2];
              α[2,1]                   (α[2,2] - ν2*k^2 - ω)]


Det_ω0(k) = det(LMat(k,0.0))

# Definition
Det_ω(k,ωr,ωi) = begin 
                  ϕ = α[1,1] + α[2,2] + (ν1 + ν2)*k^2
                  dr = Det_ω0(k) + ωr^2 - ωi^2 - ωr*ϕ
                  di = ωi*(2.0*ωr - ϕ)
                  return [dr; di]
                end  
#------------------------------ 
Det_ω(1.0,1.0,1.0)

k = 0.01
M = LMat(k,0.0)

n     = 1000

kall = range(start=0.0, stop=10.0, length = n)

# n = length(νall)

Ωr    = zeros(Float64,n,2)
Ωi    = zeros(Float64,n,2)

for i in 1:n
  local ν1 = 5.0
  local ν2 = 0.1
  local ki = kall[i]
  er,ei    = evs(ki,α,ν1,ν2)

  Ωr[i,1] =  er[1]
  Ωr[i,2] =  er[2]

  Ωi[i,1] =  ei[1]
  Ωi[i,2] =  ei[2]
end   
      

if (@isdefined h2)
  close(h2)
end
if (@isdefined h3)
  close(h3)
end

cm2         = get_cmap("RdBu_r")
h2          = figure(num=2)
pcm         = plot(kall,Ωr[:,1])
ax2         = h2.gca()
ax2.set_xlabel(L"k", fontsize=lafs)
ax2.set_ylabel(L"Ω_{r}", fontsize=lafs)

h3          = figure(num=3)
pcm         = plot(kall,Ωi[:,1])
ax3         = h3.gca()
ax3.set_xlabel(L"k", fontsize=lafs)
ax3.set_ylabel(L"Ω_{i}", fontsize=lafs)

Lkν(κ,μ1,μ2) =  [(α[1,1] - μ1*κ^2)                α[1,2];
                  α[2,1]              (α[2,2] - μ2*κ^2)]

M0 = Lkν(0.0,0.0,0.0)
t0 = tr(M0)
d0 = t0^2 - 4.0*det(M0)

println("Trace of M0: $t0")
println("Discriminant of M0: $d0")
println("Eigenvalues of M0: $(eigvals(M0))")
display(α)
println("\n")

k  = 1.0
Mk = Lkν(k,ν1,ν2)
tk = tr(Mk)
dk = tk^2 - 4.0*det(Mk)

println("Trace of Mk: $tk")
println("Discriminant of Mk: $dk")
println("Eigenvalues of Mk: $(eigvals(Mk))")







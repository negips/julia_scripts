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

#lafs = 16
#include("select_nullclines.jl")
#set               = sets[1]
#pars              = GetNullClineParams(set) 

#@variables x y
#
#Dx = Differential(x)
#Dy = Differential(y)
#fcx = zeros(Float64,4)
#fcy = zeros(Float64,3)
#fcx[1] = pars.fc0
#for i in 1:length(pars.fcx)
#  fcx[i+1] = pars.fcx[i]
#end  
#for i in 1:length(pars.fcy)
#  fcy[i] = pars.fcy[i]
#end  
#f = fcx[1] + fcx[2]*x + fcx[3]*x^2 + fcx[4]*x^3 + fcy[1]*y + fcy[2]*y^2 + fcy[3]*y^3
#
#gcx = zeros(Float64,4)
#gcy = zeros(Float64,3)
#gcx[1] = pars.gc0
#for i in 1:length(pars.gcx)
#  gcx[i+1] = pars.gcx[i]
#end  
#for i in 1:length(pars.gcy)
#  gcy[i] = pars.gcy[i]
#end  
#g = gcx[1] + gcx[2]*x + gcx[3]*x^2 + gcx[4]*x^3 + gcy[1]*y + gcy[2]*y^2 + gcy[3]*y^3
#
#SymJac = Symbolics.jacobian([g,f],[x,y]) 
#SJac   = substitute(SymJac, Dict(x => 0.0, y => 0.0), fold=false)
#Jac    = Symbolics.value.(SJac)

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

k = 1.0
M = LMat(k,0.0)

νall = range(start=0.0, stop=1.0, length = 100)
n = length(νall)

Ωr    = zeros(Float64,n,n,2)
Ωi    = zeros(Float64,n,n,2)

for i in 1:n
  for j in 1:n
     local ν1 = νall[i]
     local ν2 = νall[j]
     er,ei = evs(k,α,ν1,ν2)
#     println(er)
#     println(ei)

     Ωr[i,j,1] =  er[1]
     Ωr[i,j,2] =  er[2]

     Ωi[i,j,1] =  ei[1]
     Ωi[i,j,2] =  ei[2]
   end
end   
      
ν12d = νall*ones(Float64,1,n)
ν22d = ones(Float64,n)*νall'

if (@isdefined h2)
  close(h2)
end
if (@isdefined h3)
  close(h3)
end

cm2         = get_cmap("RdBu_r")
h2          = figure(num=2)
pcm         = pcolormesh(ν12d,ν22d,Ωr[:,:,1])
pcm.set_cmap(cm2)
ax2         = h2.gca()
colorbar()
ax2.set_xlabel(L"ν_{1}", fontsize=lafs)
ax2.set_ylabel(L"ν_{2}", fontsize=lafs)
ax2.set_title(L"Ω_{r}", fontsize=lafs)

h3          = figure(num=3)
pcm         = pcolormesh(ν12d,ν22d,Ωi[:,:,1])
pcm.set_cmap(cm2)
ax3         = h3.gca()
colorbar()
ax3.set_xlabel(L"ν_{1}", fontsize=lafs)
ax3.set_ylabel(L"ν_{2}", fontsize=lafs)
ax3.set_title(L"Ω_{i}", fontsize=lafs)

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







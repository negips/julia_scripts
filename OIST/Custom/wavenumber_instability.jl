#!/bin/julia
using LinearAlgebra
using Roots
using PyPlot

include("NullClines.jl")
include("NullClineFcn.jl")
include("ElementOf.jl")
include("NullClineParams.jl")
include("MoveFigure.jl")


include("build_nullcline_fcns.jl")
# Wavenumber Instability

function evs(k,α,ν1,ν2)
  @assert size(α) == (2,2)

  LMat =  [(α[1,1] - ν1*k^2)   α[1,2];
           α[2,1]             (α[2,2] - ν2*k^2)]
  ee   = eigvals(LMat)

  er   = real.(ee)
  ei   = imag.(ee)

  return er,ei
end  

α      = zeros(Float64,2,2)
α[1,1] = pars.gcx[1]
α[1,2] = pars.gcy[1]

α[2,1] = pars.fcx[1]
α[2,2] = pars.fcy[1]

ν1  = 0.01
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

νall = range(start=0.0, stop=10000.0, length = 1000)

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






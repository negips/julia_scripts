println("Rayleigh Oscillator. Multiple scale analysis")


import Pkg

using Blink
using PyPlot,Colors,PyCall
using LinearAlgebra
using Random

# x'' + 2*αx' + ω²x = ϵx' - ϵ/3(x')³      

ω0 = 1.0
α  = -0.5          # negative is unstable
ϵ  = 0.01
λ0 = -α + im*ω0*sqrt(1 - (α/ω0)^2)

order = 10
λ     = zeros(ComplexF64,order+1)
λ[1]  = λ0
for i in 1:order
  j = i+1
  λ[j] = λ[j-1]

  for k in 2:(j-1)
    λ[j] = λ[j] - λ[k]*λ[j-k+1]
  end
  λ[j] = λ[j]/(2.0*(λ0 + α))
end  

Ω  = 0.0 + 0.0im
λc = zeros(ComplexF64,order+1)
for i in 1:order+1
  global Ω,λc
  Ω = Ω + ϵ^(i-1)*λ[i]
  λc[i] = ϵ^(i-1)*λ[i]
end  

display(λ0)
display(Ω)




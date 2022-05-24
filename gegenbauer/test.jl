println("Generating Gegenbauer")

using LinearAlgebra
using PolynomialBases
using PyCall,PyPlot

include("SpecLib.jl")

# Generate Quadrature points
N           = 16                                # polynomial degree
lx1         = N+1                               # No of points
Basis       = LobattoLegendre(N)                # Polynomial Basis
z           = Basis.nodes

Nd          = 100                               # polynomial degree
lxd         = Nd+1                              # No of points
Basisd      = LobattoLegendre(Nd)               # Polynomial Basis
zd          = Basisd.nodes


p           = zeros(lxd)
pd          = zeros(lxd)

α = 5.0
β = 5.0

for n in 9:9
  for i in 1:lxd
    x          = zd[i]
    p[i],pd[i] = SpecLib.jacobf(n,α,β,x)
  end
  plot(zd,p)
end  




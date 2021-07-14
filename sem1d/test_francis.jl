# Testing the Francis Algorithm

include("BulgeChase.jl")
using LinearAlgebra
using Random
using PyPlot

rng = MersenneTwister(1234)

close("all")

n = 24
A = randn(rng,ComplexF64,n,n)
h = hessenberg(A)
H = convert(Matrix,h.H)

F = eigen(H)
ind = sortperm(abs.(F.values))
λ   = F.values[ind]

m   = 3

μ  = λ[1:m]
#μ  = F.values[1:m]
nμ = length(μ)

Q  = Matrix{typeof(A[1,1])}(1.0I,n,n)

Hn = copy(H)

niter = 1400
res   = zeros(Float64,niter)     # Vector
res1  = zeros(Float64,niter)     # Vector

for i in 1:niter
  global Q, Hn, res

  H0,Q0 = BulgeChase(Hn,μ,nμ)    # = Q0'*H*Q0

  Q     = Q*Q0
  
  hn    = hessenberg(H0)
  Qn    = convert(Matrix,hn.Q)
  Hn    = convert(Matrix,hn.H)

  Q     = Q*Qn

  res[i]  = abs(Hn[n-m+1,n-m])

end

#
semilogy(res)
#semilogy(res1)

Hn[n-m:n,n-m:n]

#x = polyHe1(H,μ,nμ)
# y = 0.0*x
# y[1] = -sqrt(x'*x)
# 
# u = x-y
# #u = u/u[1]
# τ = u'*x
# 
# Q = I - (1.0/τ)*u*u'

#Q0*x

# x = 0.0*x
# x[1] = im
# x[2] = 1.0
# xr   = real.(x)
# xi   = imag.(x)
# 
# 
# u = y-x
# u = u/norm(u)
# Q = I - 2.0*u*u'
# 
# Q*x



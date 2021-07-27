# Testing the Francis Algorithm

using LinearAlgebra
using Random
using PyPlot


include("BulgeChase.jl")
include("ExplicitShiftedQR.jl")
include("BulgeChase.jl")

rng = MersenneTwister(1234)

close("all")

n = 100

λs = randn(rng,ComplexF64,n)

AR = randn(rng,ComplexF64,n,n)
A  = inv(AR)*diagm(λs)*AR
h = hessenberg(A)
H = convert(Matrix,h.H)

ind = sortperm(abs.(λs),rev=false)
λ   = λs[ind]

m   = 1

μ  = λ[1:m]
#μ  = [0.]
nμ = length(μ)

Q  = Matrix{typeof(A[1,1])}(1.0I,n,n)

Hn = copy(H)

niter = 1
res   = zeros(Float64,niter)     # Vector
res1  = zeros(Float64,niter)     # Vector

#for i in 1:niter
#  global Q, Hn, res

  i = 1

#  Hpert,Q0 = CreateBulge(Hn,μ,nμ)    # = Q0'*p(H)*Q0 
  for j in 1:1
    global μ, nμ, H1, Q1
    nμ        = 20
    μ         = λ[1:nμ]
    H1,Q1     = FrancisSeq(H,μ,nμ)
  end
  display(H1[n-nμ+1,n-nμ])

  nμ        = 20
  μ         = λ[1:nμ]
  H2,Q2     = ExplicitShiftedQR(H,μ,nμ,2)
  display(H2[n-nμ+1,n-nμ])


#  μ         = λ[2:2]
#  nμ        = 1
#  H2,Q2     = FrancisAlg(H1,μ,nμ)
#  μ         = λ[3:3]
#  nμ        = 1
#  H3,Q3     = FrancisAlg(H2,μ,nμ)

#  H22,Q22   = FrancisSeq(H,λ[1:2],2)
#

  nμ        = 20
  μ         = λ[1:nμ]    
  H3,Q3     = RevFrancisSeq(H,μ,nμ)

  display(H3[nμ+1,nμ])

#  μ         = λ[1:4]
#  nμ        = 2
#  w1        = enpolyH(H,μ,nμ)
#
#  Q,w,τ     = AdjointReflectorZeros(w1',n-1,n);
#
#  A,Q0      = CreateLowerBulge(H,μ,nμ);         # A = Q0'*H*Q0








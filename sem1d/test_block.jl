#  Testing Block Orthogonal Iteration/QR/Arnoldi
println("Testing Block methods")

using LinearAlgebra
using Random

include("ExplicitShiftedQR.jl")
include("BulgeChase.jl")
include("RK4.jl")

rng = MersenneTwister(1235)

n = 300     # Matrix size
b = 3       # band size

B = randn(rng,n,n)

A = B'*B

eA = eigvals(A)

H,Q = BandedHessenberg(A,b);

eH  = eigvals(H)

#display(H)


nμ = 20 
μ  = eH[1:nμ]
ngs = 2 

H2,Q2,R2 = ExplicitShiftedQR(H,μ,nμ,ngs)
β2 = H2[n-nμ+1,n-nμ]

eH2   = eigvals(H2)

H3 = copy(H)
for j in 1:nμ
  global H3,Q3,R3
  Q3,R3 = qr(H3 - μ[j]*I)
  H3    = R3*Q3 + μ[j]*I
end  
β3 = H3[n-nμ+1,n-nμ]

println("β2:$β2")
println("β3:$β3")
println("done")

#H4    = R2*Q2 + μ[1]*I

#[eA eH eH2]

#  Testing Block Orthogonal Iteration/QR/Arnoldi
println("Testing Block methods")

using LinearAlgebra
using Random
using Printf

include("ArnIRst.jl")
include("ExplicitShiftedQR.jl")
include("BulgeChase.jl")
include("RK4.jl")

rng = MersenneTwister(1235)

VT  = Complex{Float64}

n = 12     # Matrix size
b = 2      # band size

B = randn(rng,VT,n,n)

#A = B'*B
A  = copy(B)

eA = eigvals(A)

H,Q = BandedHessenberg(A,b);

eH  = eigvals(H)

#display(H)

nμ = 1 
μ  = eH[1:nμ]
ngs = 2 

H2,Q2 = ExplicitShiftedQR(H,μ,nμ,ngs)
β2 = H2[n-nμ+1,n-nμ]
println("β2:$β2")

H3,Q3 = FrancisSeq(H,b,μ,nμ)
β3    = H3[n-nμ+1,n-nμ]
println("β3:$β3")

# Manually Doing it
Hn = copy(H)
for i in 1:nμ
  global H0,Q0,Hn,Qn
  H0,Q0 = CreateBulge(Hn,b,μ,nμ)

  Hn,Qn = ChaseBulgeDown(H0,b)
end
Q4 = (Qn*Q0)'
β4    = Hn[n-nμ+1,n-nμ]
println("β4:$β4")

# Super manual
Ht = copy(H0)
for i in 1:1
  global Qi,w,τ,Ht,k1,x
  local B
  x  = Ht[:,i]
  k1 = i+b
  Qi,w,τ = CreateReflectorZeros(x,k1,n)

  B  = Qi*Ht
  Ht = B*(Qi')
end
Qm = (Qi*Q0)'

j=1;
display([Q2[:,j] Q3[:,j] Q4[:,j] Qm[:,j]])


println("done")














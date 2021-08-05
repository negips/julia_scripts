#  Testing Block Orthogonal Iteration/QR/Arnoldi
println("Testing Block methods")

using LinearAlgebra
using Random

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

nμ = 4 
μ  = eH[1:nμ]
ngs = 2 

H2,Q2 = ExplicitShiftedQR(H,μ,nμ,ngs)
β2 = H2[n-nμ+1,n-nμ]
println("β2:$β2")

a = randn(rng,VT,n,b)
q,r = qr(a);
U   = q[:,1:b]
AV  = copy(U)

k   = 4
V   = zeros(VT,n,b*(k+1))
V[:,1:b] = copy(U)
H   = zeros(VT,b*(k+1),b*k)

for i=1:2 #k
  global V,H,AV,h

  for j in 1:b
    v  = V[:,(i-1)*b+j]
    av = A*v
    h = V[:,1:(i*b+j-1)]'*av
    H[1:(i*b+j-1),(i-1)*b+j] = h
    av .= av .- (V[:,1:(i*b+j-1)]*h)
    vnorm = norm(av)
    V[:,(i*b+j)] = av/vnorm
    H[i*b+j,(i-1)*b+j] = vnorm
  end  
  
end  




println("done")














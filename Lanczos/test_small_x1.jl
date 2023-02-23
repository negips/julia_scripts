println("Testing Fix for a small subdiagonal element")

include("NegiEig.jl")

using LinearAlgebra
using Random

n = 8
d  = randn(n)
dl = randn(n-1)
du = randn(n-1)

T = Tridiagonal(dl,d,du);
A = Matrix(T)

c = 2
r = c+2
A[r,c] = rand()

S = zeros(size(A))
S[r-1,c] = -A[r-1,c]

i = 2
B   = copy(A)
w,v = SimilarityTransform!(B,i,n)

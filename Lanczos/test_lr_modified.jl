println("Test modified LR")

using LinearAlgebra
using Random

n = 10;

el = ComplexF64
d  = randn(el,n)
dl = randn(el,n-1)
du = randn(el,n-1)

one = el(1)
zro = el(0)

M  = Matrix(Tridiagonal(dl,d,du))

#M[1,2] = zro
M[2,1] = zro

μ = eigvals(M)

en = zeros(el,n);
e1 = copy(en);
e1[1] = one'

λ = M[1];

x = (M - λ*I)*e1

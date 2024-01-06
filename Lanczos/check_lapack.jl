using LinearAlgebra
using Pseudospectra

include("MyLapack.jl")

n = 8
A = Pseudospectra.grcar(n)
F = schur(A)

T = copy(F.Schur)

T1 = copy(T)
Q = Matrix(Diagonal(ones(Float64,n)))

wantq = 1
j1    = 1
n1    = 2
n2    = 2

MyLapack.dlaexc!(wantq,j1,n1,n2,T,Q);

T2    = copy(T1)
Q2    = Matrix(Diagonal(ones(Float64,n)))

wantq = 1
j1    = 3
n1    = 2
n2    = 2
MyLapack.dlaexc!(wantq,j1,n1,n2,T2,Q2);

println("Done")

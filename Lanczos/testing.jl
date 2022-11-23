println("Testing Reductions")

using LinearAlgebra
using Random

include("HessenbergReduction.jl")
include("NegiEig.jl")
include("SetZero.jl")
include("BulgeChase.jl")

vt = Float64
tol = 1000*eps(vt)
ϵ = 1.0e-8

n  = 8
A  = randn(vt,n,n)
H  = copy(A) 
vh = UpperHessenbergReduction!(H)
SetZero!(H,tol)
T  = copy(H)
v1,w1 = UpperHessenbergtoTriDiagonal!(T);
SetZero!(T,tol)

B     = copy(T)

λ     = B[1,1] + ϵ
w0,v0 = SimilarityTransformBulge!(B,λ,1) 

C     = copy(T)
D0,Q0 = CreateBulge(C,1,λ,1)
SetZero!(D0,tol)

D     = copy(D0)

i = 1
q = ChaseBulgeDownOneStep!(D,i)
SetZero!(D,tol)

E   = copy(D)
v,w = ChaseUpperBulgeObliqueOneStep!(E,i)
SetZero!(E,tol)

i   = 2
F   = copy(F)
q   = ChaseUpperBulgeDownOneStep!(F,i)
SetZero!(F,tol)

G = copy(F)
q = ChaseBulgeDownOneStep!(G,i)
SetZero!(G,tol)

#for i in 2:0 #n-2
##  v,w = ChaseUpperBulgeObliqueOneStep!(H0,i)
##  SetZero!(H0,tol)
##  v,w = ChaseLowerBulgeObliqueOneStep!(H0,i)
##  SetZero!(H0,tol)
#
#  q1 = ChaseBulgeDownOneStep!(H0,i)
#  SetZero!(H0,tol)
#  #i = 2
#  #q = ChaseBulgeDownOneStep!(H0,i)
#  v1,w1 = ChaseUpperBulgeObliqueOneStep!(H0,i)
#  SetZero!(H0,tol)
#
#end


println("Done")






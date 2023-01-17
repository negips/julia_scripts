println("Testing Reductions")

using LinearAlgebra
using Random
using PyPlot

include("HessenbergReduction.jl")
include("NegiEig.jl")
include("SetZero.jl")
include("BulgeChase.jl")

vt = Float64
tol = 1000*eps(vt)
ϵ = 1.0e-8

n  = 12
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
v,w   = HybridBulgeChase!(D)

#i = 1
#q = ChaseBulgeDownOneStep!(D,i)
#SetZero!(D,tol)
#
#E   = copy(D)
#v,w = ChaseUpperBulgeObliqueOneStep!(E,i)
#SetZero!(E,tol)

#i   = 2
#F   = copy(E)
#q   = ChaseUpperBulgeDownOneStep!(F,i)
#SetZero!(F,tol)
#
#G = copy(F)
#q = ChaseBulgeDownOneStep!(G,i)
#SetZero!(G,tol)

#for i in 2:n-2
##  v,w = ChaseUpperBulgeObliqueOneStep!(H0,i)
##  SetZero!(H0,tol)
##  v,w = ChaseLowerBulgeObliqueOneStep!(H0,i)
##  SetZero!(H0,tol)
#
#  q1 = ChaseBulgeDownOneStep!(E,i)
#  SetZero!(E,tol)
#  #i = 2
#  #q = ChaseBulgeDownOneStep!(H0,i)
#  v1,w1 = ChaseUpperBulgeObliqueOneStep!(E,i)
#  SetZero!(E,tol)
#
#end

close("all")
h1 = figure(num=1)
s = spy(H,precision=1.0e-10,marker="o",markersize=6.0)
ax = gca()
ax.xaxis.set_visible(false)
ax.yaxis.set_visible(false)
savefig("spyfigure.png")

println("Done")






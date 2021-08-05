#  Testing Block Orthogonal Iteration
println("Testing Orthogonal Iterations")

using LinearAlgebra
using Random
using PyPlot
using PyCall

include("ExplicitShiftedQR.jl")
include("BulgeChase.jl")
include("RK4.jl")

close("all")

rng = MersenneTwister(1235)
VT  = Complex{Float64}

n = 12     # Matrix size
b = 2      # band size

B = randn(rng,VT,n,n)

#A = B'*B
A = copy(B)

eA = eigvals(A)
i1 = sortperm(real.(eA),rev=true)   # Descending order
ee = eA[i1];
eb = ee[1:b];

rcParams = PyPlot.PyDict(PyPlot.matplotlib."rcParams")
rcParams["markers.fillstyle"] = "none"
pλ0 = plot(real.(eb),imag.(eb), linestyle="none", marker="o",markersize=8)

U0 = randn(rng,VT,n,b)

# Orthogonalize
U,R = qr(U0)
Ub = U[:,1:b]
V  = copy(Ub)
V2 = zeros(Vt,n,3*b)
V2[:,1:b] = copy(V)
H  = zeros(VT,4*b,2*b)

nsteps = 10000
dt     = 0.001
ev     = zeros(VT,nsteps,b)
ev2    = zeros(VT,nsteps,b)

for i in 1:nsteps
  global V,ev
  global V2,ev2

# Simultaneous Iteration  
  for j in 1:b
    v  = V[:,j]
    v  = RK4!(A,v,dt)

#   Orthogonalize w.r.t older vectors      
    if j>1
      h = V[:,1:j-1]'*v
      v .= v .- (V[:,1:j-1]*h)
    end
    v = v./norm(v)
    V[:,j] = v
  end

# Accelerated Simultaneous Iteration  
  for j in 1:b
    v  = V2[:,j]
    v  = RK4!(A,v,dt)

#   Orthogonalize w.r.t older vectors      
    h = V2[:,1:j-1]'*v
    v .= v .- (V2[:,1:j-1]*h)
    v = v./norm(v)
    V[:,j] = v
  end


  Ar = V'*A*V
  e  = eigvals(Ar)
  ii = sortperm(real.(e),rev=true)  # Descending order
  ev[i,:] = e[ii]

end


for i in 1:b
  pλ1 = plot(real.(ev[:,i]),imag.(ev[:,i]), linestyle="-")
end  








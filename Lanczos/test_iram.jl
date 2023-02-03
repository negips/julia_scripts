# Testing IRAM implementation
println("Testing IRAM implementation")

using LinearAlgebra
using Random
using Pseudospectra
using Printf
using PyPlot

include("ArnUpd.jl")
include("ArnIRst.jl")
include("ExplicitShiftedQR.jl")
include("IRAM.jl")

close("all")

rng = MersenneTwister(1234)

vt = ComplexF64
#vt = Float64

n = 60     # Matrix size

λ  = randn(rng,vt,n)
#λ .= λ.^3
λm = diagm(0 => λ)

U = randn(rng,vt,n,n)
u = qr(U)
UQ = u.Q

#A = inv(UQ)*λm*UQ
A = Pseudospectra.grcar(n)
AT = A';

λ = eigvals(A);
λr = real.(λ)
ind = sortperm(λr,rev=true)

Nev   = 10             # Number of eigenvalues to calculate
EKryl = 10            # Additional size of Krylov space
LKryl = Nev + EKryl   # Total Size of Krylov space    

V     = zeros(vt,n,LKryl+1)
Hes   = zeros(vt,LKryl+1,LKryl)
Bg    = ones(Float64,n)             # Weight vector

v     = randn(rng,vt,n)
v     = randn(rng,vt,n)

ngs     = 2       # Number of Gram-Schmidt
nkryl   = 0
β       = 0.      # Arnoldi Residual norm
b       = 1       # Block size

h,β,r  = ArnUpd(V,b,Bg,v,nkryl,ngs)
V[:,1] = r
nkryl  = 1

ifconv = false


for i = 1:Nev
  global V,Hes,nkryl
  local h,β,r,v

  r = V[:,nkryl]
  v = A*r
  h,β,r  = ArnUpd(V,b,Bg,v,nkryl,ngs)
  Hes[1:nkryl,nkryl] = h
  Hes[nkryl+1,nkryl] = β
  V[:,nkryl+1]       = r
  nkryl              = nkryl + 1

end  


#while ~ifconv

# Major Iterations
for mi in 1:1 #20
  global V,Hes,nkryl,ifconv
  local U,G
  local β

  for i in Nev+1:LKryl

    local h,β,r,v

    r = V[:,nkryl]
    v = A*r

    V,Hes,nkryl,β,mi2 = IRAM!(V,Hes,Bg,v,nkryl,LKryl,mi,Nev,ngs)
#     nkryl = nk
#     V   = U
#     Hes = G
#     println(nkryl)

  end

#  println(nkryl)
#  U,G,nk,β,mi2 = IRAM!(V,Hes,Bg,v,nkryl,LKryl,mi,Nev,ngs)
#  nkryl = nk
#  println(nkryl)


#  U,G,nkryl,ifconv = ArnIRst(V,Hes,Bg,nkryl,LKryl+1,Nev,ngs)
#  println(["β = $β; nkryl=$nkryl"])

#  V   = U
#  Hes = G

#  println(Hes[1,1])
  β   = abs(Hes[Nev+1,Nev])
#  println("β = $β; Iteration=$mi")

  if β < 1.0e-12
    break
  end
 
end  

Ht = Hes[1:Nev,1:Nev]
F = eigen(Ht)
#println("Estimated Eigenvalues")
#display(F.values)
#
ind2 = sortperm(real.(F.values),rev=true)

display(F.values[ind2])

plot(imag.(λ),real.(λ),linestyle="none", marker="o")
plot(imag.(F.values),real.(F.values),linestyle="none",marker="s")

#λe, ϕe = eigs(A, nev=20, ncv=30, which=:LR, tol=0.0, maxiter=300)
#plot(imag.(λe),real.(λe),linestyle="none", marker="x")

#errnorm = norm(λ[ind[1:Nev]]-F.values[ind2]) 
#display("Err Norm: $errnorm")

println("Done.")




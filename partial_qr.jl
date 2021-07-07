println("Testing Partial QR")

using LinearAlgebra

ifortho = false
n  = 40; 
ev = rand(Float64,n)
evs = sort(ev,rev=true)
#evs[1] = evs[1] + 5.0


Q1 = rand(Float64,n,n)

# Orthogonalize
if ifortho
  q = copy(Q1[:,1])
  α = sqrt(q'q)
  Q1[:,1] = q/α
  for i in 2:n
     global Q1
     local q
     q = Q1[:,i]
     h = Q1[:,1:i-1]'*q
     q = q .- Q1[:,1:i-1]*h
     β = sqrt(q'q)
     Q1[:,i] = q/β
  end   
end

In = Matrix{Float64}(1.0I,n,n)
# norm(Q1'*Q1 - In)
#

A = inv(Q1)*diagm(ev)*Q1
A0 = copy(A)
p = 6

Q = Matrix{Float64}(1.0I,n,n)

niter = 300
for i in 1:niter
  global A,Q,R
  local q = A[:,1]
  local α = sqrt(q'*q)
  Q[:,1] = q/α
  R      = copy(A)
  R[:,1] = zeros(Float64,n)
  R[1,1] = α
  for j in 2:p
    local q = A[:,j];
    g = zeros(Float64,n)
    h = Q[:,1:j-1]'*q
    g[1:j-1] = h
    q = q - Q[:,1:j-1]*h
   
    h = Q[:,1:j-1]'*q
    g[1:j-1] = g[1:j-1] .+ h
    q = q - Q[:,1:j-1]*h
    β = sqrt(q'*q)
    g[j]     = β

    R[:,j]   = g
    Q[:,j]   = q/β
  end  
 
  A = R*Q


end

Q2 = Q[:,1:p]

Ar = Q2'*Q*R*Q2

F  = eigen(Ar)

ind = sortperm(real.(F.values),rev=true)

display([F.values[ind] evs[1:p]])






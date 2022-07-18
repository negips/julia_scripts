#!/bin/bash

println("Testing Jacobi Davidson Method with Harmonic Ritz values and Orthogonalization")

using LinearAlgebra
using Random
using PyPlot,Blink,PyCall
using IterativeSolvers

close("all")

n = 300            # Matrix size

Random.seed!(0)

A = randn(ComplexF64,n,n)*1.0

ei = eigvals(A)

h1 = figure(num=1,figsize=[6.0,5.0])
ax1 = subplot()
p1 = ax1.plot(real.(ei),imag.(ei),linestyle="none",marker="o")

ϵ = rand(ComplexF64)*0.1

λ   = ei[38] + ϵ         # Using this as the target eigenvalue

p1 = ax1.plot(real.(λ),imag.(λ),linestyle="none",marker="o",color="red")


u0 = randn(ComplexF64,n)

v = u0/norm(u0);

m = 30      # Search Space

V = zeros(ComplexF64,n,m)
W = zeros(ComplexF64,n,m)
H = zeros(ComplexF64,m,m)
L = zeros(ComplexF64,m,m)

w      = A*v
wnorm  = norm(w)
w      = w/wnorm
W[:,1] = w
v      = v/wnorm
V[:,1] = v
l11    = w'*v
L[1,1] = l11
h11    = w'*w
H[1,1] = h11
θ      = h11/l11
u      = v
r      = w - θ*u
t      = zeros(ComplexF64,n)

err    = 1.0
tol    = 1.0e-6
dist   = 5.0e-1

globit = 0
maxglobit = 100          # Global iterations
maxit  = 100             # GMRES iterations
while err>tol && globit<maxglobit        # Restarts

  global V,W,H,L,t,r,u,v,w,θ
  global err,globit

  globit = globit + 1

  for k in 1:m-1
    normi = 1.0/(w'*u)
    B = (I - normi*u*w')*(A - θ*I)*(I - normi*u*w')
    t = gmres(B,-r,maxiter=maxit,verbose=false)

#    Lk = L[1:k,1:k]
#    β  = inv(Lk)*W[:,1:k]'*t
#    t  = t .- V[:,1:k]*β
#    v  = t/norm(t)    
#    V[:,k+1] = v                    # v_k+1

    w             = A*t
    β             = W[:,1:k]'*w
    w1            = w - W[:,1:k]*β
    β             = W[:,1:k]'*w1
    wtilde        = w1 - W[:,1:k]*β         # Second GS
    local wnorm   = norm(wtilde)
    wtilde        = wtilde/wnorm
    W[:,k+1]      = wtilde

    Lk = L[1:k,1:k]
    β             = W[:,1:k]'*t
    ttilde        = t - V[:,1:k]*β
#    β             = W[:,1:k]'*t1              # Second GS
#    ttilde        = t1 - V[:,1:k]*β

    V[:,k+1]      = ttilde/wnorm                   # v_k+1

    lk1           = wtilde'V[:,1:k+1]              # L_k+1 row vector
    L[k+1,1:k+1]  = lk1

#   For the Orthogonal basis H should become I    
    hk1           = wtilde'*W[:,1:k+1]        # H_k+1 row vector
    H[k+1,1:k+1]  = hk1
    H[1:k+1,k+1]  = hk1'                      # H_k+1 column vector

    Lk1           = L[1:k+1,1:k+1]
    Htilde        = inv(Lk1)*H[1:k+1,1:k+1]

    f             = eigen(Htilde)
    ii            = sortperm(abs.(f.values .- λ))
    θ             = f.values[ii[1]]
    if (abs(θ-λ)>dist)
      θ = λ
    end
   
    s  = f.vectors[:,ii[1]];
    u  = V[:,1:k+1]*s
    w  = A*u
    local wnorm = norm(w)
    w  = w/wnorm
    u  = u/wnorm
    r  = w - θ*u
   
  end

  v      = u
  V[:,1] = v
  W[:,1] = w
  local l11    = w'*v
  L[1,1] = l11
  local h11    = w'*w
  H[1,1] = h11
  θ      = h11/l11
  if (abs(θ-λ)>dist)
    θ = λ
  end

  err = norm(r)
  println("Global Iteration=$globit/$maxglobit, Error Norm=$err")

end

# println("Global Iteration=$globit, Error Norm=$err")

F  = eigen(inv(L)*H) 
p2 = ax1.plot(real.(F.values),imag.(F.values),linestyle="none",marker="s",color="black")



println("Done")










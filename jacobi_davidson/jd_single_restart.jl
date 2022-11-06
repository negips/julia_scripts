#!/bin/bash

println("Testing Jacobi Davidson Method")

using LinearAlgebra
using Random
using PyPlot,PyCall
using IterativeSolvers

close("all")

n = 300            # Matrix size

Random.seed!(0)

A = randn(ComplexF64,n,n)

ei = eigvals(A)

h1 = figure(num=1,figsize=[6.0,5.0])
ax1 = subplot()
p1 = ax1.plot(real.(ei),imag.(ei),linestyle="none",marker="o")

ϵ = rand(ComplexF64)*0.2

λ   = ei[25] + ϵ         # Using this as the target eigenvalue

p1 = ax1.plot(real.(λ),imag.(λ),linestyle="none",marker="o",color="red")


u0 = randn(ComplexF64,n)

v = u0/norm(u0);

m = 10      # Search Space

V = zeros(ComplexF64,n,m)
W = zeros(ComplexF64,n,m)
H = zeros(ComplexF64,m,m)

V[:,1] = v
w      = A*v
W[:,1] = w
h11    = v'*w
H[1,1] = h11
θ      = h11
u      = v
r      = w - θ*u
t      = zeros(ComplexF64,n)

err    = 1.0
tol    = 1.0e-6

globit = 0
maxglobit = 1
while err>tol && globit<=maxglobit        # Restarts

  global V,W,H,t,r,u,v,w,θ
  global err,globit

  globit = globit + 1

  for k in 1:m-1
    B = (I - u*u')*(A - θ*I)*(I - u*u')
    t = gmres(B,-r,reltol=1.0e-6)
  
    β = V[:,1:k]'*t
    t = t .- V[:,1:k]*β
    v = t/norm(t)    
  
    V[:,k+1] = v
    w        = A*v
    W[:,k+1] = w
    hk       = V[:,1:k+1]'*w
    H[1:k+1,k+1] = hk
    hkT      = v'*W[:,1:k]
    H[k+1,1:k] = hkT
  
    f = eigen(H[1:k+1,1:k+1])
    ii = sortperm(abs.(f.values .- λ))
    θ  = f.values[ii[1]]
    if (abs(θ-λ)>0.5)
      θ = λ
    end
    s  = f.vectors[:,ii[1]];
    u  = V[:,1:k+1]*s
    w  = A*u
    r  = w - θ*u
   
  end

  v      = u
  V[:,1] = v
  W[:,1] = w
  local h11    = v'*w
  H[1,1] = h11
  θ  = h11
  if (abs(θ-λ)>0.5)
    θ = λ
  end

  err = norm(r)
  println("Global Iteration=$globit/$maxglobit, Error Norm=$err")

end

# println("Global Iteration=$globit, Error Norm=$err")

F  = eigen(H) 
p2 = ax1.plot(real.(F.values),imag.(F.values),linestyle="none",marker="s",color="black")



println("Done")










#!/bin/julia

println("Generating Plots for my Arnoldi Talk")

ENV["MPLBACKEND"]="qt5agg"

using LinearAlgebra
using Printf
using Random
using PyPlot

close("all")

ifnormal = false
lafs     = 16
lgfs     = 12

el  = ComplexF64
rng = MersenneTwister(1236) 

one = el(1)
zro = el(0)
ϵ   = rand(el)
ϵ   = one*1.0e-3
δ   = one*1.0e-4

n   = 2 
λ   = randn(rng,el,n)
λ[1] = -0.90 + 0.8im
λ[2] = λ[1] + ϵ 

A0  = Matrix(Diagonal(λ))
A   = copy(A0)

if !ifnormal
  A[1,2] = 1 - zro*ϵ
end  

B = randn(rng,el,n,n)
U, R  = qr(B)

A = U'*A*U

Fr = eigen(A)
Fl = eigen(A')

v1 = Fr.vectors[:,1]
v2 = Fr.vectors[:,2]

w1 = Fl.vectors[:,1]
w2 = Fl.vectors[:,2]


inp = v1'*v2
θ   = acos(abs(inp))

λ1  = Fr.values[1]
λ2  = Fr.values[2]

λr  = real.(λ)
λi  = imag.(λ)

@printf "Inner Product: %10.10e\n" abs(inp)
@printf "Angle θ: %10.10e\n" abs(θ)

N      = 100
λpert1 = zeros(el,N) 
λpert2 = zeros(el,N) 

for i = 1:N
  v    = v2 .+ δ*randn(el,n)
  v    = v./norm(v)

  ray1 = v'*A*v/(v'*v)
  λpert1[i] = ray1

  w    = w2 .+ δ*randn(el,n)
  wv   = w'*v
  w    = w./wv
 
  ray2 = w'*A*v/(w'*v)
  λpert2[i] = ray2

end  
#markersize=8,markeredgewidth=2,fillstyle="none"
#@printf "Rayleigh Quotient2: %4.4e\n" ray2

if ifnormal
  tit = "Normal Matrix"
  outfile = "normal.jpg"
else
  tit = "Non-Normal Matrix"
  outfile = "non_normal.jpg"
end


h1 = figure(num=1)
plot(real(λpert1),imag(λpert1),linestyle="none",marker="o",markersize=5, color="blue",markeredgewidth=2,fillstyle="none", label="Orthogonal")

#plot(real(λpert2),imag(λpert2),linestyle="none",marker="o",markersize=5, color="red",markeredgewidth=2,fillstyle="none", label="Oblique")

plot(λr,λi,linestyle="none",marker="x",markersize=10, color="black",markeredgewidth=2,fillstyle="none")

xlabel(L"λ_{r}",fontsize=lafs)
ylabel(L"λ_{i}",fontsize=lafs)
title(tit)
dλ = abs(λ2 - λ1)
xlim1 = minimum(real.(λ)) - dλ/2
xlim2 = maximum(real.(λ)) + dλ/2
#legend(fontsize=lgfs)
ax = gca();

ax.set_ylim((0.7997, 0.8003))

h1.set_figwidth(8.0)
h1.set_figheight(6.0)

@printf "Done"

savefig(outfile)





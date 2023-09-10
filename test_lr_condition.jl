println("Check Condition number of lower triangular matrices")


ENV["MPLBACKEND"]="qt5agg"

using LinearAlgebra
using PyPlot
using Random
using Printf

close("all")

n = 100

A =  [1.0 0.0;
      0.0 1.0]

condition = ones(n,n)

d = LinRange(-6.0, 0.0, n)
l = LinRange(-6.0, 0.0, n)

for i in 1:n
  for j in 1:n
    global condition
    A[2,1] = exp(l[j])
    A[1,1] = exp(d[i])

    condition[i,j] = log(cond(A))
  end
end

h1 = figure(num=1)
surf(d,l,condition,cmap="hot")
ax = gca()
#ax.set_xscale("log")
#ax.set_yscale("log")
xlabel(L"d")
ylabel(L"l")
colorbar


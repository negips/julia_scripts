println("Plot the null-clines of Cristobal's model")

using LinearAlgebra
using PyPlot,PyCall,Colors
using Polynomials

#Model:
#
# ̇q = q(-r + (δ + r)(αu^2 - (q-1)²))
#
# ̇u = -μ(1 + ρq(q - 2r))u - u³

# Null-cline equations
#
println("q null-cline: q(-r + (δ + r)(αu^2 - (q-1)²)) = 0")
println("u null-cline: -μ(1 + ρq(q - 2r))u - u³ = 0")

close("all")

lafs = 12

# Model parameters
r = 2.0           # Reynolds number
α = 1.0
δ = 1.0
ρ = 5.0
μ = -1.0

# q null-cline

nptsq       = 5000 
U           = range(-2.5, step=0.001,length=nptsq+1)
nullq       = zeros(Float64,3*(nptsq+1))
paramu      = copy(nullq)
polyq       = Polynomial([0.,0.,0.,1.],:q)
rootsq      = roots(polyq)
nrq         = 0

for i in 1:length(U)
  global nrq
  global nullq,paramu
  global polyq,rootsq
  u   = U[i]
  qc0 = 0.0
  qc1 = (δ+r)*α*u*u -r - 1.0
  qc2 = 2.0
  qc3 = -1.0
  polyq = Polynomial([qc0, qc1, qc2, qc3], :q)
  rootsq = roots(polyq)

  for j in 1:3
    if isreal(rootsq[j])
      if (abs(rootsq[j])>0.0)
        nrq = nrq + 1
        nullq[nrq]  = real(rootsq[j])
        paramu[nrq] = u
      end  
#      scatter(u,rootsq[j],color="red")
    end  
  end
end

scq = scatter(paramu[1:nrq],nullq[1:nrq],color="red")
scq.set_sizes([2])

# u null-cline
nptsu       = 10000
Q           = range(-5.0, step=0.001,length=nptsu+1)
nullu       = zeros(Float64,3*(nptsu+1))
paramq      = copy(nullu)
polyu       = Polynomial([0.,0.,0.,1.],:q)
rootsu      = roots(polyu)
nru         = 0


for i in 1:length(Q)
  global nru
  global nullu,paramq
  global polyu,rootsu
  q   = Q[i]
  uc0 = 0.0
  uc1 = -μ*(1 + ρ*q*(q - 2*r))
  uc2 = 0.0
  uc3 = -1.0
  polyu = Polynomial([uc0, uc1, uc2, uc3], :u)
  rootsu = roots(polyu)

  for j in 1:3
    if isreal(rootsu[j])
      if (abs(rootsu[j])>0.0)
        nru = nru + 1
        nullu[nru]  = real(rootsu[j])
        paramq[nru] = q
      end  
    end  
  end
end

scu = scatter(nullu[1:nru],paramq[1:nru],color="black")
scu.set_sizes([2])

ax = gca()
ax.set_ylabel(L"q",fontsize=lafs)
ax.set_xlabel(L"u",fontsize=lafs)

println("Done")
#display(rootsq)
#display(rootsu)








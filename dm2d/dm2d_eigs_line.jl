#!/usr/bin/julia

# ̇x1 = s1*x1 + 1*x2  + x1*x2
# ̇x2 = 0*x1  + s2*x2 - x1*x1

println("Pointwise Normal Eigenvalues of Dauchot-Manneville model")

using LinearAlgebra
using PyPlot

close("all")


s1 = -0.1875
s2 = -1.0

Δ  = 1.0 -4.0*s1*s2
xe = [1/2*(-1.0 + sqrt(Δ)); 1.0/(4.0*s2)*(-1.0 + sqrt(Δ))^2]
xt = [1/2*(-1.0 - sqrt(Δ)); 1.0/(4.0*s2)*(-1.0 - sqrt(Δ))^2]



Nx = 500

x1st = -0.5
x1en = 0.5

x1    = LinRange(x1st,x1en,Nx)
m     = 1.0                     # slope
x2st  = -0.7
x2    = m*(x1 .- x1st) .+ x2st  

gradf = zeros(Float64,2,2)

#gradf[1,1] = s1+x2[j]
#gradf[1,2] = 1.0+x1[i]
#gradf[2,1] = -2.0*x1[i]
#gradf[2,2] = s2
#
#println("x1=$(x1[i]); x2=$(x2[j])")
#
#ev = eigvals(gradf)
#
#display(ev)

ev1  = zeros(ComplexF64,Nx)
ev2  = zeros(ComplexF64,Nx)
evn1 = zeros(ComplexF64,Nx)
evn2 = zeros(ComplexF64,Nx)

evt1 = zeros(ComplexF64,Nx)
evt2 = zeros(ComplexF64,Nx)


for i in 1:Nx
  global ev1,ev2,evn1,evn2,evt1,evt2
  gradf[1,1] = s1+x2[i]
  gradf[1,2] = 1.0+x1[i]
  gradf[2,1] = -2.0*x1[i]
  gradf[2,2] = s2

  evals = eigvals(gradf)
  sp    = sortperm(abs.(evals))
  ev1[i] = evals[sp[2]]
  ev2[i] = evals[sp[1]]

  xdot = [s1*x1[i]+1.0*x2[i]; 0.0*x1[i]+s2*x2[i]] .+  [x1[i]*x2[i]; -x1[i]*x1[i]]
  xdot = xdot/(norm(xdot))

  Fn    = (I - xdot*(xdot'))*gradf # *(I -xdot*(xdot'))
  evals = eigvals(Fn)
  sp    = sortperm(abs.(evals))
  evn1[i] = evals[sp[2]]
  evn2[i] = evals[sp[1]]

  Ft    = xdot*(xdot')*gradf
  evals = eigvals(Ft)
  sp    = sortperm(abs.(evals))
  evt1[i] = evals[sp[2]]
  evt2[i] = evals[sp[1]]

end
#set_cmap("hsv")

# Plot the eigenvalues
h1  = figure(num=1,figsize=[8.,6.]);
pf1 = plot(x1,real.(evn1))
pf2 = plot(x1,real.(ev1))

#colorbar(ax1)

# Plot the eigenvalues
h2  = figure(num=2,figsize=[8.,6.]);
pf  = plot(x1,x2)
sc1 = scatter(xe[1],xe[2],s=20,facecolor="none",edgecolor="black")
sc2 = scatter(xt[1],xt[2],s=20,color="black")
sc3 = scatter(0.0,0.0,s=20,color="black"); draw()
ax2 = gca()
ax2.set_xlim([-1.5,2.0])
ax2.set_ylim([-1.0,1.0])

println("Done")






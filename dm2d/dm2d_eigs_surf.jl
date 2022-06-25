#!/usr/bin/julia

# ̇x1 = s1*x1 + 1*x2  + x1*x2
# ̇x2 = 0*x1  + s2*x2 - x1*x1

println("Pointwise Normal Eigenvalues of Dauchot-Manneville model")

using LinearAlgebra
using PyPlot

close("all")

Nx1 = 1000
Nx2 = 1000

x1  = LinRange(-1.5,2.0,Nx1)
x2  = LinRange(-1.0,1.0,Nx2)

s1 = -0.1875
s2 = -1.0

#i = 10
#j = 20

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

ev1  = zeros(ComplexF64,Nx1,Nx2)
ev2  = zeros(ComplexF64,Nx1,Nx2)
evn1 = zeros(ComplexF64,Nx1,Nx2)
evn2 = zeros(ComplexF64,Nx1,Nx2)

evt1 = zeros(ComplexF64,Nx1,Nx2)
evt2 = zeros(ComplexF64,Nx1,Nx2)


for i in 1:Nx1
  for j in 1:Nx2
    global ev1,ev2,evn1,evn2,evt1,evt2
    gradf[1,1] = s1+x2[j]
    gradf[1,2] = 1.0+x1[i]
    gradf[2,1] = -2.0*x1[i]
    gradf[2,2] = s2

    evals = eigvals(gradf)
    sp    = sortperm(abs.(evals))
    ev1[i,j] = evals[sp[2]]
    ev2[i,j] = evals[sp[1]]

    xdot = [s1*x1[i]+1.0*x2[j]; 0.0*x1[i]+s2*x2[j]] .+  [x1[i]*x2[j]; -x1[i]*x1[i]]
    xdot = xdot/(norm(xdot))

    Fn    = (I - xdot*(xdot'))*gradf
    evals = eigvals(Fn)
    sp    = sortperm(abs.(evals))
    evn1[i,j] = evals[sp[2]]
    evn2[i,j] = evals[sp[1]]

    Ft    = xdot*(xdot')*gradf
    evals = eigvals(Ft)
    sp    = sortperm(abs.(evals))
    evt1[i,j] = evals[sp[2]]
    evt2[i,j] = evals[sp[1]]

  end
end

Δ  = 1.0 -4.0*s1*s2
xe = [1/2*(-1.0 + sqrt(Δ)); 1.0/(4.0*s2)*(-1.0 + sqrt(Δ))^2]
xt = [1/2*(-1.0 - sqrt(Δ)); 1.0/(4.0*s2)*(-1.0 - sqrt(Δ))^2]

set_cmap("hsv")

# Plot the eigenvalues
h1  = figure(num=1,figsize=[8.,6.]);
sf  = surf(x1,x2,real.(ev1),cmap=get_cmap(),linewidth=0.0)
ax1 = gca(); ax1.view_init(90,-90); draw()
sc1 = scatter3D(xe[1],xe[2],[2.0],s=20,facecolor="none",edgecolor="black")
sc2 = scatter3D(xt[1],xt[2],[2.0],s=20,color="black")
sc3 = scatter3D(0.0,0.0,[2.0],s=20,color="black"); draw()
#colorbar(ax1)








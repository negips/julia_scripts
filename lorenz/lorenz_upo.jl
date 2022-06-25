#!/usr/bin/julia

println("UPOs of the Lorenz system")

# ̇x = σ(y - x)
# ̇y = ρx -y -zx
# ̇z = -βz + xy

# Standard parameters:
# σ = 10.0
# ρ = 28.0        # also R
# β = 8.0/3.0

using LinearAlgebra
using Random
using PyPlot,Blink


include("Lorenz_ddt.jl")
include("LinLorenz_ddt.jl")

include("Lorenz_nstep.jl")
include("LinLorenz_nstep.jl")

include("RK4_lorenz.jl")
include("RK4_linlorenz.jl")


close("all")

lafs = 16

σ = 10.0
ρ = 28.0          # also R
β = 8.0/3.0

T  = 10.0
nsteps = 100000
dt = T/nsteps

nflds = 4               # no of intermediate fields

nstep_seg = Int(nsteps/nflds)

newton_fac = 0.001        # do only a partial newton step
niter   = 10000             # No of Iterations
plotupd = 10

Random.seed!(12)        # initialize seed
X0    = rand(Float64,3,nflds).*10.0
X     = zeros(Float64,3,nflds)
dx    = zeros(Float64,3)
res   = zeros(Float64,3)
rnorm = zeros(Float64,niter)

x = abs(X0[1,1])
y = abs(X0[2,1])
z = abs(X0[3,1])
nn = 20000
h1  = figure(num=1,figsize=[8.,6.])
ax1 = subplot()
cm  = get_cmap("tab20")

x,y,z = Lorenz_nstep!(x,y,z,σ,ρ,β,dt,nn)

nn = nstep_seg
for i =1:nflds
  global X,x,y,z
  x,y,z = Lorenz_nstep!(x,y,z,σ,ρ,β,dt,nn)
  X[1,i] = x
  X[2,i] = y
  X[3,i] = z

  p1 = ax1.scatter(x,z,color=cm(i-1),s=20,marker="s")
end


display(X)
scarr = []
plarr = []

for i in 1:niter
   
  global X,dx,res,rnorm,ax1,scarr

  if mod(i,plotupd)==0
    global scarr,plarr
    for sc in scarr
      sc.remove()
    end
    for pl in plarr
      pl[1].remove()
    end  
  end  


  for j in 1:nflds

#   Starting point      
    x0 = X[1,j]
    y0 = X[2,j]
    z0 = X[3,j]
    if (mod(i,plotupd)==0 || i==1)
      sc = ax1.scatter(x0,z0,color=cm(j-1),s=20)
      if i==1
        push!(scarr,sc)
      else
        scarr[j] = sc
      end
    end  

#   Target point      
    l = mod(j,nflds)+1
    xt = X[1,l]
    yt = X[2,l]
    zt = X[3,l]

#   Starting point    
    xi = X[1,j]
    yi = X[2,j]
    zi = X[3,j]

    xhis = zeros(Float64,nstep_seg)
    zhis = zeros(Float64,nstep_seg)
    for is = 1:nstep_seg
      xi,yi,zi = RK4_lorenz!(xi,yi,zi,σ,ρ,β,dt)
      xhis[is] = xi
      zhis[is] = zi
    end
    if (mod(i,plotupd)==0 || i==1)
      plhis = ax1.plot(xhis,zhis,color=cm(j-1))
      if i==1
        push!(plarr,plhis)
      else
        plarr[j] = plhis
      end  
    end 

#   Residuals
    res = zeros(Float64,3)   
    res[1] = xi-xt
    res[2] = yi-yt
    res[3] = zi-zt

    rnorm[i] = norm(res)

    v   = rand(Float64,3)
    nor = norm(v)
    v   = v/nor

    V       = zeros(Float64,3,3)
    Hess    = zeros(3,3)
    V[:,1]  = v

    for k = 1:3
      x0i = X[1,j] 
      y0i = X[2,j]
      z0i = X[3,j] 
      xp = copy(v[1])
      yp = copy(v[2])
      zp = copy(v[3])
      xp,yp,zp,x0i,y0i,z0i = LinLorenz_nstep!(xp,yp,zp,x0i,y0i,z0i,σ,ρ,β,dt,nstep_seg)
      v[1]   = xp
      v[2]   = yp
      v[3]   = zp
      if k==1
        α      = V[:,1]'*v
        Hess[1,k] = α
        v      = v - V[:,1]*α
        nor    = norm(v)
        Hess[k+1,k] = nor
        v      = v/nor
        V[:,k+1] = v
      elseif k==2
        α2     = V[:,1:k]'*v
        Hess[1:k,k] = α2
        v      = v - V[:,1:k]*α2
        nor    = norm(v)
        Hess[k+1,k] = nor
        v      = v/nor
        V[:,k+1] = v
      elseif k==3
        α3     = V[:,1:k]'*v
        Hess[1:k,k] = α3
        v      = v - V[:,1:k]*α3
        nor    = norm(v)      # this should be zero now
      end
    
    end           # Krylov subspace building

    if (mod(i,plotupd)==0)
      pause(0.001)
    end  

#   Newton Step
    dx    = -inv(Hess)*V'*res
    X[:,j] = X[:,j] .+ dx*newton_fac

  end       # nflds    

end         # i = 1:niter 

#plot(rnorm)

# plot(x,z)
# ax1 = gca();
# ax1.set_xlabel(L"x",fontsize=lafs)
# ax1.set_ylabel(L"z",fontsize=lafs)

display(X)

println("Done.")







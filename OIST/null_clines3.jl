#!/bin/julia

println("Calculate the null clines for Eqns 2.1a, 2.1b in H. Meinhardt (1995) The algorithmic beauty of sea shells")

using LinearAlgebra
using Random
using Printf
using PyPlot,PyCall
using Roots

include("Meinhardt.jl")
include("GetBDF.jl")
include("GetEXT.jl")

close("all")

lafs = 16

ra = 1.0
s  = 1.0
ba = 0.1

ca = ra/(2*s)

rb = 1.0
bb = 1.0
cb = rb/s

s  = 1

Adotxy(x,y) = Meinhardt_25(x,y,s)[1]
Adotx(x) = Adotxy(x,b)
Adoty(y) = Adotxy(a,y)

Bdotxy(x,y) = Meinhardt_25(x,y,s)[2]
Bdotx(x) = Bdotxy(x,b)
Bdoty(y) = Bdotxy(a,y)


# Initial point
arange0 = -40.0
arange1 = +40.0
b0 = 25.000
b  = b0
aa = find_zeros(Adotx,arange0,arange1)
a  = aa[1]

nsteps = 50000

nroots   = length(aa)

aroots_a = zeros(nsteps,nroots)
aroots_b = zeros(nsteps,nroots)
ϕ        = zeros(nsteps,nroots)

dτ = -1.0e-3

h1 = figure(num=1)

# Find Activator null-cline
for n in 1:nroots
  global a,b,a0,b0,da,db

  b   = b0
  global an  = find_zeros(Adotx,arange0,arange1);
  b   = b0 + dτ
  global an1 = find_zeros(Adotx,arange0,arange1);

  if n==1
    da = an1[1] - an[1]
    a  = an[1]
    b  = b0
  else
    da = an1[2] - an[2]
    a  = an[2]
    b  = b0
  end  
  db = dτ

  for i in 1:nsteps

    θa = atan(dτ,da)
    θb = atan(dτ,db)
    θ      = atan(db,da)  
    ϕ[i,n] = θ
  
    if abs(da) > abs(db)
      a  = a + da
      b1 = find_zero(Adoty,b,atol=1.0e-8);
      db = b1 - b
      b  = b1
    else
      b  = b + db
      a1 = find_zero(Adotx,a,atol=1.0e-08);
      da = a1 - a
      a  = a1
    end
  
    aroots_a[i,n] = a
    aroots_b[i,n] = b
  end
end  

plot(aroots_b,aroots_a,linestyle="-")

#

# Initial point
arange0 = -40.0
arange1 = +40.0
b0 = -40.000
b  = b0
aa = find_zeros(Bdotx,arange0,arange1)
a  = aa[1]

nroots = length(aa)

broots_a = zeros(nsteps,nroots)
broots_b = zeros(nsteps,nroots)

dτ = 1.0e-3

# Find Inhibitor null-cline
for n in 1:nroots
  global a,b,a0,b0,da,db

  b   = b0
  global an  = find_zeros(Bdotx,arange0,arange1);
  b   = b0 + dτ
  global an1 = find_zeros(Bdotx,arange0,arange1);

  if n==1
    da = an1[1] - an[1]
    a  = an[1]
    b  = b0
  else
    da = an1[2] - an[2]
    a  = an[2]
    b  = b0
  end  
  db = dτ

  for i in 1:nsteps

    θa = atan(dτ,da)
    θb = atan(dτ,db)
    θ      = atan(db,da)  
    ϕ[i,n] = θ
  
    if abs(da) > abs(db)
      a  = a + da
      b1 = find_zero(Bdoty,b,atol=1.0e-8);
      db = b1 - b
      b  = b1
    else
      b  = b + db
      a1 = find_zero(Bdotx,a,atol=1.0e-08);
      da = a1 - a
      a  = a1
    end
  
    broots_a[i,n] = a
    broots_b[i,n] = b
  end
end  

plot(broots_b,broots_a,linestyle="--")

ax1 = h1.axes[1]
ax1.set_xlabel(L"b", fontsize=lafs)
ax1.set_ylabel(L"a", fontsize=lafs)
#ax1.legend(fontsize=12)















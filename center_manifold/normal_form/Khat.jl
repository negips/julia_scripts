#!/bin/julia

# Build Khat

println("Building Khat")

using FileIO
using PyPlot

DEV   = load("direct_GL_5.jld2");
V     = get(DEV,"evec",[])
v     = V[:,5];
xg    = get(DEV,"xg",[])
Bg    = get(DEV,"Bg",[])

AEV   = load("adjoint_GL_5.jld2");
W     = get(AEV,"evec",[])
w     = W[:,5];

v1    = [v;zero(v)] 
v2    = [zero(v); conj.(v)]

w1    = [w;zero(w)]
w2    = [zero(w); conj.(w)]

n     = 2

close("all")

cm    = get_cmap("tab10")
plot(xg,real.(v),linestyle="-",color=cm(0))
plot(xg,imag.(v),linestyle="--",color=cm(0))

plot(xg,real.(w),linestyle="-",color=cm(1))
plot(xg,imag.(conj.(w)),linestyle="--",color=cm(1))

x0    = 7.0
ξ     = exp.(-(xg .- x0).^2)
β     = sqrt(ξ'*(Bg.*ξ))
ξ     = ξ./β

plot(xg,ξ,linestyle="-",color=cm(2),linewidth=2)

# Forcings
h           = 5
f1          = (1.0/sqrt(2.0))*[ξ;ξ]
f2          = [ξ;       zero(ξ)]
f3          = [zero(ξ); ξ]
f4          = [ξ;       zero(ξ)]
f5          = [zero(ξ); ξ]
Bg2         = [Bg;Bg]

ΓH          = zeros(ComplexF64,n,h)

ΓH[1,2]     = w1'*(Bg2.*f2)

ΓH[2,3]     = w2'*(Bg2.*f3)

display(ΓH[1,2])
display(ΓH[2,3])

println("Done.")






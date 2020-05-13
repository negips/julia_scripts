println("Eigen spectrum of convection diffusion operator")

using LinearAlgebra
using PyPlot,PyCall

close("all")

# Include the function files
include("sem_main.jl")

r,c = size(L);
Lsub = L[2:r,2:c];

F  = eigen(Lsub);

F.values

lr = real.(F.values);
li = imag.(F.values);

plot(li,lr,linestyle="none",marker="o",markersize=8)

#u = sin.(xm1[:,1]);
#plot(xm1[:,1],u)
#
#du = gradx[:,:,1]*u;
#plot(xm1[:,1],du)
#
#xd = intpm1d*xm1[:,1];
#ud = intpm1d*u;
#plot(xd,ud)




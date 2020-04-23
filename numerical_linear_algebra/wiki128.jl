using Pkg
#Pkg.add("PyPlot")
#Pkg.add("LaTeXStrings")
#Pkg.add("PyCall")

using LinearAlgebra,MAT # for: eig, norm, etc
using MatrixDepot, Random,TimerOutputs,PyPlot,LaTeXStrings
using PyCall

#include("arnoldi.jl")
#include("GS_npass.jl")

println("Wiki, Problem 1-36")

lmbda = [10;3;1];

PyCall.PyDict(matplotlib["rcParams"])["markers.fillstyle"] = "full"
PyCall.PyDict(matplotlib["rcParams"])["axes.labelsize"] = 16;
PyCall.PyDict(matplotlib["rcParams"])["axes.grid"] = true;
PyCall.PyDict(matplotlib["rcParams"])["grid.linestyle"] = ":";

pl = plot(real(lmbda),imag(lmbda),linestyle="none", marker="o",markersize=10)
ylabel(L"\lambda_i")
xlabel(L"\lambda_r")

PyCall.PyDict(matplotlib["rcParams"])["markers.fillstyle"] = "none"
pl = plot(real(lmbda[1]),imag(lmbda[2]),linestyle="none", marker="o",markersize=24)

c = (lmbda[2]+lmbda[3])/2;
cx = real(c);
cy = imag(c);
rad = abs(lmbda[2]-lmbda[3])/2;

theta = [i for i in 0:0.001:2*pi];
x = cx .+ rad.*cos.(theta);
y = cy .+ rad.*sin.(theta);

dist=abs(lmbda[1]-c);

conv = abs(rad/dist);

#figure()
pl2 = plot(x,y, linestyle="--",color="k");
pl3 = plot(cx,cy,linestyle="none",color="k",linewidth=4,marker="+",markersize=12)

savefig("1-29.png")



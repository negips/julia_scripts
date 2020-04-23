using Pkg
#Pkg.add("PyPlot")
#Pkg.add("LaTeXStrings")
#Pkg.add("PyCall")

using LinearAlgebra,MAT # for: eig, norm, etc
using MatrixDepot, Random,TimerOutputs,PyPlot,LaTeXStrings
using PyCall

include("arnoldi.jl")
include("GS_npass.jl")

println("Homework 1, Exercise 6a")

B = matread("Bwedge.mat")["B"];

lmbda = eigvals(B);

PyCall.PyDict(matplotlib["rcParams"])["markers.fillstyle"] = "none"
PyCall.PyDict(matplotlib["rcParams"])["axes.labelsize"] = 16;
PyCall.PyDict(matplotlib["rcParams"])["axes.grid"] = true;
PyCall.PyDict(matplotlib["rcParams"])["grid.linestyle"] = ":";

pl = plot(real(lmbda),imag(lmbda), linestyle="none", marker=".");
ylabel(L"\lambda_{i}");
xlabel(L"\lambda_{r}");

# Making an approximate circle
im_lambda = imag(lmbda);

indmax = argmax(im_lambda);
cx = real(lmbda[indmax])-5;
cy = 0;
rad = imag(lmbda[indmax])+2;

theta = [i for i in 0:0.001:2*pi];
x = cx .+ rad.*cos.(theta);
y = cy .+ rad.*sin.(theta);

#figure()
pl2 = plot(x,y, linestyle="--",color="k");
pl3 = plot(cx,cy,linestyle="none",color="k",linewidth=4,marker="+",markersize=12)




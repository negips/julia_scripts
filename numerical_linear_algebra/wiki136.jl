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

lmbda = [3+2im;-2+2im;-1-4im];
mu = -1+im;

shift = lmbda .-mu;
shift_invert = 1 ./shift;

PyCall.PyDict(matplotlib["rcParams"])["markers.fillstyle"] = "full"
PyCall.PyDict(matplotlib["rcParams"])["axes.labelsize"] = 16;
PyCall.PyDict(matplotlib["rcParams"])["axes.grid"] = true;
PyCall.PyDict(matplotlib["rcParams"])["grid.linestyle"] = ":";

pl = plot(real(shift_invert),imag(shift_invert),linestyle="none", marker="o",markersize=10)
ylabel(L"\lambda_i(A-\mu I)^{-1}")
xlabel(L"\lambda_r(A-\mu I)^{-1}")

PyCall.PyDict(matplotlib["rcParams"])["markers.fillstyle"] = "none"
pl = plot(real(shift_invert[2]),imag(shift_invert[2]),linestyle="none", marker="o",markersize=24)

c = (shift_invert[1]+shift_invert[3])/2;
cx = real(c);
cy = imag(c);
rad = abs(shift_invert[3]-shift_invert[1])/2;

theta = [i for i in 0:0.001:2*pi];
x = cx .+ rad.*cos.(theta);
y = cy .+ rad.*sin.(theta);

dist=abs(shift_invert[2]-c);

conv = abs(rad/dist);

#figure()
pl2 = plot(x,y, linestyle="--",color="k");
pl3 = plot(cx,cy,linestyle="none",color="k",linewidth=4,marker="+",markersize=12)

savefig("shift_invert.png")



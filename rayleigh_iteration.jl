##
#
println("Rayleigh Iteration for Eigenvalue problem")

using LinearAlgebra
#using Gadfly
using PyPlot,PyCall


close("all")

# v0 = [1. 2. 2.5];
# n  = length(v0);

#n  = 10;
#v0 = rand(Float64,n);
#A  = rand(Float64,n,n);

vi = [0.01 0.1 0.2]*0.1;
vr = [1.3 1.1 -1.15]*0.1;

n = length(vr);
n2 = 2*n;

i = sqrt(Complex(-1.));

#A = zeros(Complex,n2,n2);
#for j in 1:n
#  global A    
#  k = 2*(j-1)+1;
#  e = exp(vr[j] + i*vi[j]);
#  A[k,k] = e;
#  k = 2*(j-1)+2;
#  e = exp(vr[j] - i*vi[j]);
#  A[k,k] = e;
#end  


A = rand(Float64,n2,n2);

m = 1;
u = rand(Float64,n2);
u = u/norm(u);

niters = 100;
eye    = Matrix(Complex(1.0)I,n2,n2);

muhis = zeros(Complex,niters);

for j in 1:niters
  global u, mu    
  mu = u'*A*u;

  muhis[j] = mu;

  u  = (A - mu.*eye)\u;
  u  = u/norm(u);
 
end  

println(mu)

plot(real.(muhis))

F = eigen(A);

dist = sqrt.(real.(F.values).^2 + imag.(F.values).^2)

Fs = eigen(0.5*(A + A'));











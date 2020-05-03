#
println("Testing Shifted simultaneous iterations");

using Random
using LinearAlgebra
using PyPlot,PyCall

close("all")

Random.seed!(111);        # Initialize seed for consistent results

vi = [0.5 1.1 1.2]*0.01;
vr = [0.1 0.2 -0.15]*0.01;

n = length(vr);

A = zeros(Complex,2*n,2*n);

i = sqrt(Complex(-1.));

for j in 1:n
  global A    
  k = 2*(j-1)+1;
  e = exp(vr[j] + i*vi[j]);
  A[k,k] = e;
  k = 2*(j-1)+2;
  e = exp(vr[j] - i*vi[j]);
  A[k,k] = e;
end  

# Rotate A
V1 = Complex.(rand(Float64,n2,n2));
V2 = i.*Complex.(rand(Float64,n2,n2));

V =  V1 + V2;

V[:,1] = V[:,1]/norm(V[:,1]);

for j in 2:n2
  global V    
  u = copy(V[:,j]);
  α = V[:,1:j-1]'*u;
  u = u - V[:,1:j-1]*α;
  u = u/norm(u);
  V[:,j] = u;
end  

A = V'*A*V;

F = eigen(A);

ee = F.values;
er = real.(ee);
ei = imag.(ee);

dist = sqrt.(er.^2 + ei.^2);

plot(ei,er,linestyle="",marker=".")

n2 = 2*n;
m  = 2;
U = Complex.(rand(Float64,n2,m));

U[:,1] = U[:,1]/norm(U[:,1]);

for j in 2:m
  global U    
  u = copy(U[:,j]);
  α = U[:,1:j-1]'*u;
  u = u - U[:,1:j-1]*α;
  u = u/norm(u);
  U[:,j] = u;
end  


# Power iterations
sigma = 0.0;
niter = 5000;

h2    = figure(num=2);

for j in 1:niter
  global A,U,Ar
  global sigma

  A1 = copy(A);
  A1 = A1 - sigma.*I;

  U = A1*U;
  U[:,1] = U[:,1]/norm(U[:,1]);
  
  for j in 2:m
    global U    
    u = copy(U[:,j]);
    α = U[:,1:j-1]'*u;
    u = u - U[:,1:j-1]*α;
    u = u/norm(u);
    U[:,j] = u;
  end  

  Ar = (i.*U)'*A*(i.*U);
  F2 = eigen(Ar);

#  plot(j,real.(F2.values[1]),linestyle="", marker=".",color="blue")
#  plot(j,real.(F2.values[2]),linestyle="", marker=".",color="blue")

  plot(j,imag.(F2.values[1]),linestyle="", marker=".",color="blue")
  plot(j,imag.(F2.values[2]),linestyle="", marker=".",color="blue")


end













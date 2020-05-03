##
#
println("Testing Rayleigh Quotients of projected Matrices")

using LinearAlgebra
using Gadfly


# v0 = [1. 2. 2.5];
# n  = length(v0);

n  = 10;
v0 = rand(Float64,n);
A  = zeros(Float64,n,n);

for i in 1:n
   global A
   A[i,i] = v0[i];
end

m = 4;
U = rand(Float64,n,m);
U[:,1] = U[:,1]/norm(U[:,1]);

# Orthogonalize
for i in 2:m
   global U   
   u  = copy(U[:,i]);
   α  = U[:,1:i-1]'*u;
   u  = u - U[:,1:i-1]*α;
   u  = u/norm(u);
   U[:,i] = u;
end 

norm(1.0I - U'*U);

Ar = U'*A*U

F = eigen(Ar);
println(F.values)


ray = U[:,1]'*A*U[:,1];
for i in 2:m
  global ray    
  r2  = U[:,i]'*A*U[:,i];
  ray = vcat(ray,r2);
end  

println(ray)




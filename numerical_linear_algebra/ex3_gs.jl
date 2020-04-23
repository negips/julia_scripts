using LinearAlgebra # for: eig, norm, etc

#import Pkg
#Pkg.add("MatrixDepot")
using Random

include("arnoldi.jl")

function GS_npass(Q,w,k,ngs)
#  write your function

#  Inputs
#  Q - Matrix of orthogonal k orthogonal vectors
#  k - no of orthogonal vectors
#  w - vector to orthogonalize
#  ngs - No of Gram-Schmidt loops

#  Outputs
#  h - components along Qk vectors
#  beta - component orthogonal to Qk
#  qortho - vector orthogonal to Qk

   rw,cl=size(Q);
   h = zeros(Float64,k,1);  # Initialize to 0

   for igs in 1:ngs
#     global h
#     global w 

     g = Q[:,1:k]'*w;
     h = h .+ g; 
     w = w .- Q[:,1:k]*g;
   end   

   wnorm = norm(w);
   qnew = w./wnorm;
   beta = w'*qnew;
   qortho = beta.*qnew;
      
#  Check orthogonality of qnew
#   test = norm(Q[:,1:k]'*qnew);
#   println("Orthogonality of new Vector: $test")

   return h,beta,qortho
end

Random.seed!(0)
A=randn(100,100); b=randn(100);
m=5;
Q,H=arnoldi(A,b,m);
should_be_zero1=norm(Q*H-A*Q[:,1:m]);
should_be_zero2=norm(Q'*Q-I);

println("Arnoldi Factorization condition satisfied to : $should_be_zero1")
println("Q Matrix Orthogonality satisfied to : $should_be_zero2")



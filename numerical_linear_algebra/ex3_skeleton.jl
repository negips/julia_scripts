using LinearAlgebra # for: eig, norm, etc
include("arnoldi.jl")

function hw1_good_gs(Q,w,k)
    # write your function
end


A=randn(100,100); b=randn(100);
m=10;
Q,H=arnoldi(A,b,m);
should_be_zero1=norm(Q*H-A*Q[:,1:m])
should_be_zero2=norm(Q'*Q-I)

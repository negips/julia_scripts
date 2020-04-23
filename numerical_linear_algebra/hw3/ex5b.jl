#!/usr/bin/julia
# Source code for Homework 3

println("Homework 3, Exercise 5b")

#import Pkg
#Pkg.add("PyPlot")
#Pkg.add("PyCall")
#Pkg.add("Blink")

using Blink,TimerOutputs
using PyPlot,PyCall,LinearAlgebra

include("schur_parlett.jl")

close("all")      # close all open figures

const to = TimerOutput();

m=100;
A = rand(m,m);
A = A/norm(A);

fsine=x->sin(x);
fadd=x->x^2
#fpower=(x,y)->x^y
# F=schur_parlett(A,fadd)

Ns = [i for i in 20:20:1400];
np = length(Ns);

timer_naive = zeros(np,1);
timer_sp = zeros(np,1);


for i in 1:np
# Naive Method

  global B
  global timer_naive
  global timer_sp
  N=Ns[i];

  B = Matrix{Float64}(I,size(A));
  B = A;
  for j=1:N-1
    B= @timeit to "Naive $N" A*B;
  end
  timer_naive[i]=TimerOutputs.time(to["Naive $N"])*1.0e-6;  # mili seconds

  fpn=(x)->x^N;
  F= @timeit to "schur_parlett $N" schur_parlett(A,fpn)
  timer_sp[i]=TimerOutputs.time(to["schur_parlett $N"])*1.0e-6;  # mili seconds

#  diff=norm(B-F);
#  println("Norm(B-F)=$diff for N=$N")    
  
end  

pl = plot(Ns[2:np],timer_naive[2:np], marker=".", label="Naive");
pl2 = plot(Ns[2:np],timer_sp[2:np], marker=".", label="Schur-Parlett")
ylabel("CPU-time (ms)",fontsize="14")
xlabel(L"N",fontsize="14")
legend()

fname="ex5b_sp.png";
savefig(fname);

to

#!/usr/bin/julia
# Source code for Homework 3

println("Homework 3, Exercise 3c")

#import Pkg
#Pkg.add("PyPlot")
#Pkg.add("PyCall")
#Pkg.add("Blink")

using Blink,TimerOutputs
using PyPlot,PyCall,LinearAlgebra

include("alpha_example.jl")
include("hessenberg_alg2.jl")

close("all")      # close all open figures

const to = TimerOutput();

ifmod = false;    # if we use modified GS
ngs = 1;          # No of GS passes

errtol = 1.0e-10; # error tolerance for sub diagonal elements

alpha = 1;
m_array = [5 10 100 200 300 400]*1;

for m in m_array
# Create Matrix A
  global A
  global H
  A=alpha_example(alpha,m);
  H = @timeit to "hessenberg $m" hessenberg_alg2(A);   

end

to

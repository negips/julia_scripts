#!/usr/bin/julia
# Source code for Homework 3

println("Homework 3, Exercise 5a")

#import Pkg
#Pkg.add("PyPlot")
#Pkg.add("PyCall")
#Pkg.add("Blink")

using Blink,TimerOutputs
using PyPlot,PyCall,LinearAlgebra

include("schur_parlett.jl")

close("all")      # close all open figures

#const to = TimerOutput();

A = [1. 4. 4.; 3. -1. 3.; -1. 4. 4.];

fsine=x->sin(x);
F=schur_parlett(A,fsine)

# to

#!/bin/julia
# Tests for Nek Fourier implementation.
println("Tests for Nek2.5D implementation.")
println("Convection Diffusion spectra")

using LinearAlgebra
using PyCall
using PyPlot


α     = 1.0
β     = 1.0
γ     = 0.750       # Spanwise

U     = 1.0
V     = 1.0
W     = 1.2

ν     = 2000;

Conv  = zeros(ComplexF64,3,3)
Diff  = zeros(ComplexF64,3,3)


Conv[1,1] = im*(U*α + V*β + W*γ)
Conv[2,2] = im*(U*α + V*β + W*γ)
Conv[3,3] = im*(U*α + V*β + W*γ)


Diff[1,1]   = 2*α^2 + β^2 + γ^2
Diff[2,2]   = α^2 + 2*β^2 + γ^2
Diff[3,3]   = α^2 + β^2 + 2*γ^2

Diff[1,2]   = α*β
Diff[1,3]   = α*γ

Diff[2,1]   = β*α 
Diff[2,3]   = β*γ

Diff[3,1]   = γ*α 
Diff[3,2]   = γ*β

Diff        = -1.0/ν*Diff

A           = -Conv .+ Diff;

λ           = eigvals(A)

nsteps      = 150.0
dt          = 0.001
B           = exp(A*dt*nsteps)

λ2          = eigvals(B)

[λ λ2]

#!/bin/julia

"""
    BarkleyPipe(q,u,r)

Compute the Following Equation:

 ∂q/∂t = q×(r + u - U0 - (r + δ)×(q - 1.0)² )

 ∂u/∂t = ϵ1×(U0 - u) + ϵ2×(Ub - u)×q

# Examples
```julia-repl
julia> adot,bdot = BarkleyPipe(1.0,1.0,0.5)
(-0.5,0.1)
```
"""
function BarkleyPipe(u,q,r)

  U0  = 2.0  # Centre line velocity
  Ub  = 1.0  # Bulk velocity

#  ζ   = 0.8  # Convection parameter
  δ   = 0.1
  ϵ   = 0.1
  κ   = 2.0
  ϵ1  = ϵ
  ϵ2  = κ*ϵ

  g   = ϵ1*(U0 .- u) .+ ϵ2*(Ub .- u).*q
  f   = q.*(r .+ u .- U0 .- (r + δ)*(q .- 1.0).^2 ) #.+ σ*q

  return [g f]
end
#---------------------------------------------------------------------- 

"""
    LiesWithin(x,x0,x1)

Check if x ∈ (x0,x1)

# Examples
```julia-repl
julia> LiesWithin(1.0,0.0,2.0)
true
```
"""
function LiesWithin(x,x0,x1)

  f = x>x0 && x<x1

  return f
end

#---------------------------------------------------------------------- 

function BarkleyNullClines(r)

  g(x,y)     = BarkleyPipe(x,y,r)[1]
  f(x,y)     = BarkleyPipe(x,y,r)[2]

  nsteps = 100000

  qr0  = -20.0
  qr1  = +20.0
  u0   =   4.0
  q0   =  u0
  τq   = -1.0e-4

  qnc_u,qnc_q = NullClines(f,u0,qr0,qr1,nsteps,τq)

  nsteps = 50000

  qr0  = -10.0
  qr1  = +10.0
  u0   =   4.0
  q0   =  u0
  τu   = -1.0e-4

  unc_u,unc_q = NullClines(g,u0,qr0,qr1,nsteps,τu)

  return qnc_q,qnc_u,unc_q,unc_u
end
#---------------------------------------------------------------------- 














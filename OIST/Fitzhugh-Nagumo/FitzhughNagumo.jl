#!/bin/julia

"""
      FitzhughNagumo(u,v,a,b)

      Compute the Following Equation:

      ∂u/∂t = 3×u - u^3 - v

      ∂v/∂t = u - a - b×v

"""
function FitzhughNagumo(u,v,a,b)


  f   = 3*u .- u.^3 .- v
  g   = u .- a .- b*v

  return [f g]
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

function FitzhughNagumoNullClines(a,b)

  f(x,y)     = FitzhughNagumo(x,y,a,b)[1]
  g(x,y)     = FitzhughNagumo(x,y,a,b)[2]

  nsteps = 100000

  vr0  = -20.0
  vr1  = +20.0
  u0   =   3.0
  v0   =  u0
  τu   = -1.0e-4

  unc_u,unc_v = NullClines(f,u0,vr0,vr1,nsteps,τu)

  nsteps = 50000

  vr0  = -10.0
  vr1  = +10.0
  u0   =  -3.5
  v0   =  u0
  τv   =  1.0e-3

  vnc_u,vnc_v = NullClines(g,u0,vr0,vr1,nsteps,τv)

  return unc_u,unc_v,vnc_u,vnc_v
end
#---------------------------------------------------------------------- 














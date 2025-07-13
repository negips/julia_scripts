#!/bin/julia

function FHN2_G(v,u,R,m)

  g   = R*(u .- m*v)

  return g
end
#---------------------------------------------------------------------- 
function FHN2_F(v,u,a1,b1)

  α1  = 0.0
  α2  = 1.0
  α3  = a1
  α4  = b1

  f   = -(u .- α1).*(u .- α2).*(u .- α3) .+ α4*v

  return f
end

#---------------------------------------------------------------------- 

"""
      FitzhughNagumo2(v,u,a,b,R,m)

      Compute the Following Equation:


      g = ∂v/∂t = u - c*v

      f = ∂u/∂t = -(u - 0.0)*(u - 1.0)*(u - a) - b*v

      Returns [g f]
"""
function FitzhughNagumo2(v,u,a1,b1,R,m)

  α1  = 0.0
  α2  = 1.0
  α3  = a1
  α4  = b1

  g   = FHN2_G(v,u,R,m)
  f   = FHN2_F(v,u,a1,b1)

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

function FitzhughNagumoNullClines2(a,b,R,m)

  g(x,y)     = FitzhughNagumo2(x,y,a,b,R,m)[1]
  f(x,y)     = FitzhughNagumo2(x,y,a,b,R,m)[2]

  nsteps = 300000

  vr0  = -20.0
  vr1  = +20.0
  u0   =  10.0
  v0   =  u0
  τu   = -1.0e-4

  unc_u,unc_v = NullClines(f,u0,vr0,vr1,nsteps,τu)

  nsteps = 50000

  vr0  = -10.0
  vr1  = +10.0
  u0   =  -6.0
  v0   =  u0
  τv   =  1.0e-3

  vnc_u,vnc_v = NullClines(g,u0,vr0,vr1,nsteps,τv)

  return vnc_u,vnc_v,unc_u,unc_v
end
#---------------------------------------------------------------------- 
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














#!/bin/julia
"""
    ElementOf(x,x0,x1)

Check if x ∈ (x0,x1)


# Examples
```julia-repl
julia> ElementOf(1.0,0.0,2.0)
true
```
"""
function ElementOf(x,x0,x1)

  f = x>x0 && x<x1

  return f
end




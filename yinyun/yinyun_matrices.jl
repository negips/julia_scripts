#!/bin/julia

# println("Testing Yinyun maximization")

function YinyunMatrix(nr::Integer,nc::Integer,α::Real,β::Real,kai::Real,δ::Real,γ::Real,Ca::Real)

  n   = nr*nc
  M   = zeros(Float64,n,n)

  b   = 0.25

  k         = [3.0*α;   2.0*α;      α;          0]
  kappa     = [0.0;     β;          2.0*b*β;    3*(b^2)*β]
  f         = [2*kai;   kai;        0.0]
  g         = [0.0;     δ;          2*b*δ]

  XI        = CartesianIndices((1:nr,1:nc))

  for ci in XI

    println(ci)

    r0 = ci[1]
    c0 = ci[2]

    r_1 = r0 - 1
    r1  = r0 + 1
    c_1 = c0 - 1
    c1  = c0 + 1

    l0   = r0  + (c0-1)*nr
    lr_1 = r_1 + (c0-1)*nr
    lr1  = r1  + (c0-1)*nr
    lc_1 = r0  + (c_1-1)*nr
    lc1  = r0  + (c1-1)*nr

    println([l0; lr_1; lr1; lc_1; lc1])

    # Across Rows
    if (r0 == 1)
      M[l0,l0]    = M[l0,l0] - g[r0]           # Going Up
      M[l0,l0]    = M[l0,l0] - Ca*f[r0]        # Going Down

      # M[l0,lr_1]  = M[l0,lr_1] + Ca*f[r_1]     # From Top
      M[l0,lr1]   = M[l0,lr1] + g[r1]          # From Bottom
    elseif (r0 == nr)
      M[l0,l0]    = M[l0,l0] - g[r0]           # Going Up
      M[l0,l0]    = M[l0,l0] - Ca*f[r0]        # Going Down

      M[l0,lr_1]  = M[l0,lr_1] + Ca*f[r_1]     # From Top
      # M[l0,lr1]   = M[l0,lr1] + g[r1]          # From Bottom
    else
      M[l0,l0]    = M[l0,l0] - g[r0]           # Going Up
      M[l0,l0]    = M[l0,l0] - Ca*f[r0]        # Going Down

      M[l0,lr_1]  = M[l0,lr_1] + Ca*f[r_1]     # From Top
      M[l0,lr1]   = M[l0,lr1] + g[r1]          # From Bottom
    end

    # Across Columns
    if (c0 == 1)
      M[l0,l0]    = M[l0,l0] - kappa[c0]       # Going Left
      M[l0,l0]    = M[l0,l0] - Ca*k[c0]        # Going Right

      # M[l0,lc_1]  = M[l0,lc_1] + Ca*k[c_1]     # From Left
      M[l0,lc1]   = M[l0,lc1] + kappa[c1]      # From Right
    elseif (c0 == nc)
      M[l0,l0]    = M[l0,l0] - kappa[c0]       # Going Left
      M[l0,l0]    = M[l0,l0] - Ca*k[c0]        # Going Right

      M[l0,lc_1]  = M[l0,lc_1] + Ca*k[c_1]     # From Left
      # M[l0,lc1]   = M[l0,lc1] + kappa[c1]      # From Right
    else
      M[l0,l0]    = M[l0,l0] - kappa[c0]       # Going Left
      M[l0,l0]    = M[l0,l0] - Ca*k[c0]        # Going Right

      M[l0,lc_1]  = M[l0,lc_1] + Ca*k[c_1]     # From Left
      M[l0,lc1]   = M[l0,lc1] + kappa[c1]      # From Right
    end

  end       # ci in XI

  return M
end  

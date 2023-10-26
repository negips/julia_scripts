#!/bin/julia
"""
    Meinhardt_21(A,B,S)

Compute the Following Equation:

 ∂A/∂t = S×A² - ra×A + ba×S

 ∂B/∂t = S×A² - rb×B + bb 

# Examples
```julia-repl
julia> adot,bdot = Meinhardt_21(2.0,1.0,1.0)
(3.0,4.0)
```
"""
function Meinhardt_21(A,B,S)

# A  - Activator
# B  - Inhibitor
# S  - Substrate

  ra = 1.0        # Linear Removal of activator
  ba = 1.0        # Constant activator production/offset term
  
  rb = 1.0        # Linear removal of inhibitor
  bb = 1.0        # Constant inhibitor production/offset term

  Adot  = S.*A.^2 ./B .- ra*A .+ ba*S
  Bdot  = S.*A.^2 .- rb*B + bb

  return Adot,Bdot
end

#---------------------------------------------------------------------- 
"""
    Meinhardt_24(A,B,S)

Compute the Following Equation:

 ϕ = A²/(1 + sa×A²) + ba

 ∂A/∂t = S×B×ϕ - ra*×A

 ∂B/∂t = bb - S×B×ϕ - rb×B

# Examples
```julia-repl
julia> adot,bdot = Meinhardt_24(3.0,1.0,1.0)
(-1.1,-1.9)
```
"""
function Meinhardt_24(A,B,S)

# A  - Activator
# B  - Inhibitor
# S  - Substrate

  sa = 1.0        # Saturation for activator
  ra = 1.0        # Linear Removal of activator
  ba = 1.0        # Constant activator production/offset term
  
  rb = 1.0        # Linear removal of inhibitor
  bb = 1.0        # Constant inhibitor production/offset term

  Astar2 = (A.^2)./(1 + sa*A.^2) .+ ba

  Adot  = S.*B.*Astar2 .- ra*A
  Bdot  = bb .- S.*B.*Astar2 .- rb*B

  return Adot,Bdot
end

#---------------------------------------------------------------------- 

function Meinhardt_25(A,B,S)

# A  - Activator
# B  - Inhibitor
# S  - Substrate

  sa = 1.0        # Saturation for activator
  ra = 1.0        # Linear Removal of activator
  ba = 1.0        # Constant activator production/offset term
  
  rb = 1.0        # Linear removal of inhibitor
  bb = 1.0        # Constant inhibitor production/offset term

  Adot  = S.*(A.^2) .- ra*B.*A .+ S*ba
  Bdot  = bb .- S.*(A.^2) .- rb*B

  return Adot,Bdot
end

#---------------------------------------------------------------------- 

function Meinhardt_26(A,B,C,S)

# A  - Activator
# B  - Inhibitor
# C  - Switching system along with A
# S  - Substrate

  sa = 1.0        # Saturation for activator
  ra = 1.0        # Linear Removal of activator
  ba = 1.0        # Constant activator production/offset term
  
  rb = 1.0        # Linear removal of inhibitor
  bb = 1.0        # Constant inhibitor production/offset term

  sc = 1.0        # Saturation constant
  rc = 1.0        # Linear removal rate for C

  Adot  = S./(sa .+ C.^2) .- ra*A .+ ba
  Bdot  = rb*A .- rb*B
  Cdot  = S./(sc .+ (A./B).^2) .- rc*C

  return Adot,Bdot,Cdot
end

#---------------------------------------------------------------------- 

function Meinhardt_31(A,B,S)

# A  - Activator
# B  - Inhibitor
# S  - Substrate

  sa = 1.0        # Saturation for activator
  ra = 1.0        # Linear Removal of activator
  ba = 1.0        # Constant activator production/offset term

  sb = 1.0        # Saturation for activator via the inhibitor in the denominator
  rb = 1.0        # Linear removal of inhibitor
  bb = 1.0        # Constant inhibitor production/offset term

  Adot  = S.*A.^2 ./((sb .+ B).*(1 .+ sa*A.^2)) .- ra*A .+ S.*ba
  Bdot  = S.*A.^2 .- rb*B + bb

  return Adot,Bdot
end

#---------------------------------------------------------------------- 

function Meinhardt_51(A,B,C,S)

# A  - Activator
# B  - Inhibitor
# C  - Switching system along with A
# S  - Substrate

  sa = 1.0        # Saturation for activator
  ra = 1.0        # Linear Removal of activator
  ba = 1.0        # Constant activator production/offset term
  
  rb = 1.0        # Linear removal of inhibitor
  bb = 1.0        # Constant inhibitor production/offset term

  sc = 1.0        # Saturation constant
  rc = 1.0        # Linear removal rate for C

  Astar2 = (A.^2 .+ ba)./(1 .+ sa*(A.^2))

  Adot  = S.*B.*Astar2./(sb .+ sc*C) .- ra*A
  Bdot  = bb .-  S.*B.*Astar2./(sb .+ sc*C) .- rb*B
  Cdot  = rc*(A .- C)

  return Adot,Bdot,Cdot
end

#---------------------------------------------------------------------- 

function Meinhardt_52(A,B,C,S)

# A  - Activator
# B  - Inhibitor
# C  - Switching system along with A
# S  - Substrate

  sa = 1.0        # Saturation for activator
  ra = 1.0        # Linear Removal of activator
  ba = 1.0        # Constant activator production/offset term
  
  rb = 1.0        # Linear removal of inhibitor
  bb = 1.0        # Constant inhibitor production/offset term
  bc = 1.0        # Constant C production

  sc = 1.0        # Saturation constant
  rc = 1.0        # Linear removal rate for C

  Astar2 = (A.^2 .+ ba)./(1 .+ sa*(A.^2))

  Adot  = S.*Astar2.*(B .+ C) .- ra*A
  Bdot  = bb .-  S.*B.*Astar2 .- rb*B
  Cdot  = bc .-  S.*C.*Astar2 .- rc*C

  return Adot,Bdot,Cdot
end

#---------------------------------------------------------------------- 

function Meinhardt_53(A,B,C,S)

# A  - Activator
# B  - Inhibitor
# C  - Switching system along with A
# S  - Substrate

  sa = 1.0        # Saturation for activator
  ra = 1.0        # Linear Removal of activator
  ba = 1.0        # Constant activator production/offset term
  
  rb = 1.0        # Linear removal of inhibitor
  bb = 1.0        # Constant inhibitor production/offset term
  bc = 1.0        # Constant C production

  sc = 1.0        # Saturation constant
  rc = 1.0        # Linear removal rate for C

  Adot  = (S./C).*((A.^2)./B .+ ba) .- ra*A
  Bdot  = bb .- rb*B .+ rb*(A.^2)./C
  Cdot  = rc*(A .- C)

  return Adot,Bdot,Cdot
end

#---------------------------------------------------------------------- 

"""
    Meinhardt_1987_61(A,B,C,D,S)

Equation 1a,b from Meinhardt & Klinger (1987) A Model for Pattern Formation on the Shells of Molluscs,
 J. theor. Biol. 126, 63-89

 A     - Activator
 B     - Inhibitor

 ϕ     = A²/(1 + κ×A²)

 ∂A/∂t = ρ×(ϕ + ρ0)/B - μ×A

 ∂B/∂t = ρ1 + ρ×ϕ - ν×B

# Examples
```julia-repl
julia> adot,bdot = Meinhardt_61(2.0,1.0)
(3.0,4.0)
```
"""
function Meinhardt_61(A,B,C,D,S)

# A  - Activator
# B  - Inhibitor
# C  - Global control of rb
# S  - Substrate
# D  - Activator saturation

  sb = 1.0        # Saturation for activator
  sd = 1.0        # Saturation for activator based on D
  ra = 1.0        # Linear Removal of activator
  ba = 1.0        # Constant activator production/offset term
  
  rb = 1.0        # Linear removal of inhibitor
  bb = 1.0        # Constant inhibitor production/offset term
  bc = 1.0        # Constant C production

#  sc = 1.0        # Saturation constant
#  rc = 1.0        # Linear removal rate for C

  rbeff = rb/C      

  Adot  = S.*((A.^2)./(sb .+ sd*D.^2 .+ B) .+ ba) .- ra*A
  Bdot  = bb .+ s*(A.^2) .- (rbeff)*B 

#  Ddot  = rd*(A .- D)

  return Adot,Bdot
end

#---------------------------------------------------------------------- 

"""
    Meinhardt_1987_1(A,S)

Equation 1a,b from Meinhardt & Klinger (1987) A Model for Pattern Formation on the Shells of Molluscs,
 J. theor. Biol. 126, 63-89

 A     - Activator
 S     - Substrate

 ϕ     = A²/(1 + κ×A²) + ρ0

 ∂A/∂t = ρ×S×ϕ - μ×A

 ∂S/∂t = σ - ρ×S×ϕ - ν×S

# Examples
```julia-repl
julia> adot,bdot = Meinhardt_1987_1(2.0,1.0,1.0)
(3.0,4.0)
```
"""
function Meinhardt_1987_1(A,S)

# A  - Activator
# S  - Substrate

  κ  = 0.0         # Saturation constant for production
#  κ = 0.08 # (Figure 1c)                  
  ρ0 = 0.001       # production offset

  μ  = 0.01        # Linear Removal of activator
  ρ  = 0.01        # Production term coefficient

  
  ν  = 0.0         # Linear removal of inhibitor
  σ  = 0.015       # Constant substrate production/offset term

  Da = 0.002       # Diffusion coefficient for A
  Ds = 0.4         # Diffusion coefficient for S

  ϕ     = (A.^2)./(1 .+ κ*A.^2) .+ ρ0

  Adot  = ρ*S.*ϕ .- μ*A
  Sdot  = σ .- ρ*S.*ϕ .- ν*S

  return Adot,Sdot
end

#---------------------------------------------------------------------- 

"""
    Meinhardt_1987_2(A,B)

Equation 1a,b from Meinhardt & Klinger (1987) A Model for Pattern Formation on the Shells of Molluscs,
 J. theor. Biol. 126, 63-89

 A     - Activator
 B     - Inhibitor

 ϕ     = A²/(1 + κ×A²)

 ∂A/∂t = ρ×(ϕ + ρ0)/B - μ×A

 ∂B/∂t = ρ1 + ρ×ϕ - ν×B

# Examples
```julia-repl
julia> adot,bdot = Meinhardt_1987_2(2.0,1.0)
(3.0,4.0)
```
"""
function Meinhardt_1987_2(A,B)

# A  - Activator
# B  - Inhibitor

  κ  = 0.15        # Saturation constant for production
#  κ  = 0.08 # (Figure 2e)
#  κ  = 0.4  # (Figure 2f)
  ρ0 = 0.01        # production offset

  μ  = 0.02        # Linear Removal of activator
  ρ  = 0.2         # Production term coefficient

  
  ν  = 0.02        # Linear removal of inhibitor
  ρ1 = 0.01        # Constant substrate production/offset term (Not specified in paper)

  Da = 0.01        # Diffusion coefficient for A
  Db = 0.4         # Diffusion coefficient for S

  ϕ     = (A.^2)./(1 .+ κ*A.^2)

  Adot  = ρ*(ϕ .+ ρ0)./B .- μ*A
  Bdot  = ρ1 .+ ρ*ϕ .- ν*B

  return Adot,Bdot
end

#---------------------------------------------------------------------- 

"""
    Meinhardt_1987_2_branching(A,B,R)

Equation 1a,b from Meinhardt & Klinger (1987) A Model for Pattern Formation on the Shells of Molluscs,
 J. theor. Biol. 126, 63-89

 A     - Activator
 B     - Inhibitor

 ϕ     = A²/(1 + κ×A²)

 ∂A/∂t = ρ×(ϕ + ρ0)/B - μ×A

 ∂B/∂t = ρ1 + ρ×ϕ - ν×B

 Parameters for Figure 8 in paper.
 κ  = 0.25        # Saturation constant for production
 ρ0 = 0.01        # production offset

 μ  = 0.02        # Linear Removal of activator
 ρ  = 0.1         # Production term coefficient

 
 ν  = 0.0014      # Linear removal of inhibitor
 ρ1 = 0.01        # Constant substrate production/offset term (Not specified in paper)

 Da = 0.015       # Diffusion coefficient for A
 Db = 0.0         # Diffusion coefficient for S


# Examples
```julia-repl
julia> adot,bdot = Meinhardt_1987_2_branching(2.0,1.0)
(3.0,4.0)
```
"""
function Meinhardt_1987_2_branching(A,B,R)

# A  - Activator
# B  - Inhibitor

  κ  = 0.25        # Saturation constant for production
  ρ0 = 0.1         # production offset

  μ  = 0.1         # Linear Removal of activator
  ρ  = 0.1         # Production term coefficient
  
  ν  = 0.0014      # Linear removal of inhibitor
  ρ1 = 0.001        # Constant substrate production/offset term (Not specified in paper)

  Da = 0.015       # Diffusion coefficient for A
  Db = 0.0         # Diffusion coefficient for S

  B0 = 0.1         # Saturation of Activator via Inhibitor

  ϕ     = (A.^2)./(1.0 .+ κ*A.^2)

  Adot  = ρ*(ϕ .+ ρ0)./(B .+ B0) .- μ*A
 
  Bdot  = ρ1 .+ ρ*ϕ .- (ν/R)*B

  return Adot,Bdot
end

#---------------------------------------------------------------------- 








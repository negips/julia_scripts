#!/bin/julia

# Sets: 1         : Pulses
#     : 2         : Slugs. No change in instability threshold
#     : 3         : Pulses. Smaller instability threshold
#     : 4         : Slugs. Smaller instability threshold
#     : 5         : Extreme Slugs. Smaller instability threshold
#     : 6         : Extreme Pulse collapse
#     : 7         : Extreme slugs
#     : 8         : Weak slugs
#     : 10        : Limit-cycline Oscillation. Activation dominated
#     : 11        : Limit-cycline Oscillation. De-activation dominated
#     : 12        : Limit-cycline Oscillation. De-activation dominated (cubic in x)
#     : 13        : Two fixed points - upper and lower branch.
#     : 14        : Symmetric fixed points - upper and lower branch
#     : 15        : Symmetric LCO
#     : 16        : Two Asymmetric fixed points
#     : 17        : Unstable G-null-cline
#     : 18        : Marginal Instability
#     : 19        : Stable Oscillatory Linear mode
#     : 51        : Dynamic Switching (λ)
#     : 52        : Dynamic Switching (λ)
#     : 53        : Dynamic Switching (λ)

sets  = [3]
println("Nullcline Set(s): $sets")


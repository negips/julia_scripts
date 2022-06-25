#!/usr/bin/julia
function Lorenz_nstep!(x,y,z,σ,ρ,β,dt,nsteps)

  for i in 1:nsteps
     
     x,y,z = RK4_lorenz!(x,y,z,σ,ρ,β,dt)
  
  end   
  
  return x,y,z

end





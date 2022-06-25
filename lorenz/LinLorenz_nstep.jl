#!/usr/bin/julia
function LinLorenz_nstep!(xp,yp,zp,x0,y0,z0,σ,ρ,β,dt,nsteps)
# Linearized Lorenz nstep  

  for i in 1:nsteps
     
     xp,yp,zp = RK4_linlorenz!(xp,yp,zp,x0,y0,z0,σ,ρ,β,dt)
     x0,y0,z0 = RK4_lorenz!(x0,y0,z0,σ,ρ,β,dt)
 
  end   
  
  return xp,yp,zp,x0,y0,z0

end





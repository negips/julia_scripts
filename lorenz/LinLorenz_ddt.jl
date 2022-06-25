# Time Derivative of the Lorenz system
function LinLorenz_ddt(xp,yp,zp,x0,y0,z0,σ,ρ,β)

  dxdt  = -σ*xp         + σ*yp            + 0.0*zp
  dydt  = (ρ-z0)*xp     - 1.0*yp          - x0*zp
  dzdt  = y0*xp         + x0*yp          - β*zp

  return dxdt,dydt,dzdt
end  


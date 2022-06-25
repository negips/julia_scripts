# Time Derivative of the Lorenz system
function Lorenz_ddt(x,y,z,σ,ρ,β)

  dxdt  = σ*(y - x)
  dydt = (ρ*x - y - z*x)
  dzdt  = (x*y - β*z)

  return dxdt,dydt,dzdt
end  


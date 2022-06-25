# 4th Order Runge-Kutta Steps for the Lorenz system
function RK4_lorenz!(x,y,z,σ,ρ,β,dt)

 dxdt,dydt,dzdt = Lorenz_ddt(x,y,z,σ,ρ,β)

# dxdt  = σ*(y - x)
# dydt  = (ρ*x - y - z*x)
# dzdt  = (β*z + x*y)

  x1 = x + dt/2.0*dxdt
  y1 = y + dt/2.0*dydt
  z1 = z + dt/2.0*dzdt

  dxdt1,dydt1,dzdt1 = Lorenz_ddt(x1,y1,z1,σ,ρ,β)

# dxdt1  = σ*(y1 - x1)
# dydt1  = (ρ*x1 - y1 - z1*x1)
# dzdt1  = (β*z1 + x1*y1)

  x2 = x + dt/2.0*dxdt1
  y2 = y + dt/2.0*dydt1
  z2 = z + dt/2.0*dzdt1

  dxdt2,dydt2,dzdt2 = Lorenz_ddt(x2,y2,z2,σ,ρ,β)

# dxdt2  = σ*(y2 - x2)
# dydt2  = (ρ*x2 - y2 - z2*x2)
# dzdt2  = (β*z2 + x2*y2)

  x3 = x + dt*dxdt2
  y3 = y + dt*dydt2
  z3 = z + dt*dzdt2

  dxdt3,dydt3,dzdt3 = Lorenz_ddt(x3,y3,z3,σ,ρ,β)

# dxdt3  = σ*(y3 - x3)
# dydt3  = (ρ*x3 - y3 - z3*x3)
# dzdt3  = (β*z3 + x3*y3)

  x  = x + dt/6.0* (dxdt + 2.0*dxdt1 + 2.0*dxdt2 +dxdt3)
  y  = y + dt/6.0* (dydt + 2.0*dydt1 + 2.0*dydt2 +dydt3)
  z  = z + dt/6.0* (dzdt + 2.0*dzdt1 + 2.0*dzdt2 +dzdt3)

  return x,y,z 
end  


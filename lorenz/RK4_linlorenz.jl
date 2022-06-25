# 4th Order Runge-Kutta Steps for the Lorenz system
function RK4_linlorenz!(xp,yp,zp,x0,y0,z0,σ,ρ,β,dt)

  dxdt,dydt,dzdt = LinLorenz_ddt(xp,yp,zp,x0,y0,z0,σ,ρ,β)

  x1 = xp + dt/2.0*dxdt
  y1 = yp + dt/2.0*dydt
  z1 = zp + dt/2.0*dzdt

  dxdt1,dydt1,dzdt1 = LinLorenz_ddt(x1,y1,z1,x0,y0,z0,σ,ρ,β)

  x2 = xp + dt/2.0*dxdt1
  y2 = yp + dt/2.0*dydt1
  z2 = zp + dt/2.0*dzdt1

  dxdt2,dydt2,dzdt2 = LinLorenz_ddt(x2,y2,z2,x0,y0,z0,σ,ρ,β)

  x3 = xp + dt*dxdt2
  y3 = yp + dt*dydt2
  z3 = zp + dt*dzdt2

  dxdt3,dydt3,dzdt3 = LinLorenz_ddt(x3,y3,z3,x0,y0,z0,σ,ρ,β)

  xp  = xp + dt/6.0* (dxdt + 2.0*dxdt1 + 2.0*dxdt2 +dxdt3)
  yp  = yp + dt/6.0* (dydt + 2.0*dydt1 + 2.0*dydt2 +dydt3)
  zp  = zp + dt/6.0* (dzdt + 2.0*dzdt1 + 2.0*dzdt2 +dzdt3)

  return xp,yp,zp
end  


# 4th Order Runge-Kutta Steps
function RK4!(λdot,a,λ,dt)

  localprec = eltype(λ)
  two       = localprec(2)
  six       = localprec(6)

  λ1 = λ + dt/two*λdot(a,λ)
  λ2 = λ + dt/two*λdot(a,λ1)
  λ3 = λ + dt*λdot(a,λ2)
  λ  = λ + dt/six*(λdot(a,λ) + two*λdot(a,λ1) + two*λdot(a,λ2) + λdot(a,λ3))

  return λ 
end  


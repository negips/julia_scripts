# 4th Order Runge-Kutta Steps
function RK4!(M,v,dt)

  localprec = eltype(v[1])
  two = localprec(2)
  six = localprec(6)

  v1 = v .+ dt/two*M*v
  v2 = v .+ dt/two*M*v1
  v3 = v .+ dt*M*v2
  v .= v .+ dt/six*(M*(v .+ two*v1 .+ two*v2 .+ v3))

  return v    
end  


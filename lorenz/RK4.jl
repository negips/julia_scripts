# 4th Order Runge-Kutta Steps
function RK4!(M,v,dt)

  v1 = v .+ dt/2.0*M*v
  v2 = v .+ dt/2.0*M*v1
  v3 = v .+ dt*M*v2
  v .= v .+ dt/6.0*(M*(v .+ 2.0*v1 .+ 2.0*v2 .+ v3))

  return v    
end  


# 4th Order Runge-Kutta Steps
# OP is a function acting on v
function OP_RK4!(OP,v,dt)

  localprec = eltype(v[1])
  two = localprec(2)
  six = localprec(6)

  v1 = v .+ dt/two*OP(v)
  v2 = v .+ dt/two*OP(v1)
  v3 = v .+ dt*OP(v2)
  v .= v .+ dt/six*(OP(v) .+ two*OP(v1) .+ two*OP(v2) .+ OP(v3))

  return v    
end  


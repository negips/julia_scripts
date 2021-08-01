# 4th Order Runge-Kutta Steps
function SEM_RK4!(v,dt,nel,OP,B,lx1,Q,QT,prec)

  vi = copy(v)
  v1 = copy(v)
  v2 = copy(v)
  v3 = copy(v)

  if prec==BigFloat
    one = BigFloat(1.0)
    two = BigFloat(2.0)
    six = BigFloat(6.0)
  else
    one = 1.0
    two = 2.0
    six = 6.0
  end  
 
  Bdssum = Q*(QT*B)
  Binv   = one./Bdssum

  for i in 1:nel
    j1 = (i-1)*lx1 
    j2 = i*lx1
    vi[j1:j2] .= OP[:,:,i]*v[j1:j2]
  end
  v1 .= v .+ (dt/two)*Binv.*Q*(QT*vi)

  for i in 1:nel
    j1 = (i-1)*lx1 
    j2 = i*lx1
    vi[j1:j2] .= OP[:,:,i]*v1[j1:j2]
  end
  v2 .= v .+ (dt/two)*Binv.*Q*(QT*vi)

  for i in 1:nel
    j1 = (i-1)*lx1 
    j2 = i*lx1
    vi[j1:j2] .= OP[:,:,i]*v2[j1:j2]
  end
  v3 .= v .+ (dt)*Binv.*Q*(QT*vi)

  for i in 1:nel
    j1 = (i-1)*lx1 
    j2 = i*lx1
    vi[j1:j2] .= OP[:,:,i]*(v[j1:j2] .+ two*v1[j1:j2] .+ two*v2[j1:j2] .+ v3[j1:j2])
  end
  v .= v .+ (dt/six)*Binv.*Q*(QT*vi)

  return v    
end  


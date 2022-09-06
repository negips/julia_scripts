# 4th Order Runge-Kutta Steps
function SEM_RK4!(v,dt,nel,lx1,OP,B,Binv,Q,QT,prec)

  vi = copy(v)
  v1 = copy(v)
  v2 = copy(v)
  v3 = copy(v)

  zro = prec(0.0)
  one = prec(1.0)
  two = prec(2.0)
  six = prec(6.0)
 
#  Bdssum = Q*(QT*B)
#  Binv   = one./Bdssum

  vi  = SEM_Local_Apply(v,nel,lx1,OP,Binv,Q,QT,prec)
  v1 .= v .+ (dt/two)*vi
  v1[1] = zro + zro*im 

  vi  = SEM_Local_Apply(v1,nel,lx1,OP,Binv,Q,QT,prec)
  v2 .= v .+ (dt/two)*vi
  v2[1] = zro + zro*im 

  vi  = SEM_Local_Apply(v2,nel,lx1,OP,Binv,Q,QT,prec)
  v3 .= v .+ (dt)*vi
  v3[1] = zro + zro*im 

#  for i in 1:nel
#    j1 = (i-1)*lx1+1
#    j2 = i*lx1
#    vi[j1:j2] .= OP[:,:,i]*(v[j1:j2] .+ two*v1[j1:j2] .+ two*v2[j1:j2] .+ v3[j1:j2])
#  end

  vi  = SEM_Local_Apply(v,nel,lx1,OP,Binv,Q,QT,prec)
  v  .= v .+ (dt/six)*vi
  vi  = SEM_Local_Apply(v1,nel,lx1,OP,Binv,Q,QT,prec)
  v  .= v .+ (dt/six)*two*vi
  vi  = SEM_Local_Apply(v2,nel,lx1,OP,Binv,Q,QT,prec)
  v  .= v .+ (dt/six)*two*vi
  vi  = SEM_Local_Apply(v3,nel,lx1,OP,Binv,Q,QT,prec)
  v  .= v .+ (dt/six)*vi

  v[1] = zro + zro*im 

  return v    
end  
#---------------------------------------------------------------------- 
function SEM_Local_Apply(v,nel,lx1,OP,Binv,Q,QT,prec)

  vi   = copy(v)
  vout = copy(v)

  if prec==BigFloat
    zro = BigFloat(0.0)
    one = BigFloat(1.0)
  else
    zro = 0.0
    one = 1.0
  end  

  for i in 1:nel
    j1 = (i-1)*lx1+1 
    j2 = i*lx1
    vi[j1:j2] .= OP[:,:,i]*v[j1:j2]
  end
  vout .= Binv.*Q*(QT*vi)
  vout[1] = zro + zro*im 

  return vout

end  




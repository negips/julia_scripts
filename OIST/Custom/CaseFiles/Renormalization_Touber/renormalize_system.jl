# function renormalize_system!(pars)
# 
#   ndec    = 5
# 
#   f(x,y)  = FXY(x,y,pars.fc0,pars.fcx,pars.fcy)
# 
#   f00(y)    = f(0.0,y)
#   roots_f   = find_zeros(f00,-10.0,10.0)
#   sort!(roots_f)
#   α1_2      = roots_f[2]      # Make second root == 1.0
# 
#   # f-parameters
# 
#   pars.fc0        = pars.fc0/α1_2
#   pars.fcx[1]     = pars.fcx[1]
# 
#   pars.fcy[1]     = pars.fcy[1]
#   pars.fcy[2]     = pars.fcy[2]*α1_2
#   pars.fcy[3]     = pars.fcy[3]*(α1_2)^2
# 
#   pars.gc0        = pars.gc0/α1_2
# 
#   return α1_2
# end 
# #---------------------------------------------------------------------- 
function renormalize_system!(pars)

  ndec    = 5

  f(x,y)  = FXY(x,y,pars.fc0,pars.fcx,pars.fcy)

  f00(y)    = f(0.0,y)
  roots_f   = find_zeros(f00,-10.0,10.0)
  sort!(roots_f)
  α1_2      = roots_f[2]      # Make second root == 1.0

  # A'      = A/Anorm
  Anorm     = α1_2
  
  # B'      = B/Bnorm
  Bnorm     = Anorm/1.0 

  # f-parameters

  pars.fc0        = pars.fc0/Anorm
  pars.fcx[1]     = pars.fcx[1]*Bnorm/Anorm

  pars.fcy[1]     = pars.fcy[1]
  pars.fcy[2]     = pars.fcy[2]*Anorm
  pars.fcy[3]     = pars.fcy[3]*(Anorm)^2

  pars.gc0        = pars.gc0/Bnorm
  pars.gcx[1]     = pars.gcx[1]
  pars.gcy[1]     = pars.gcy[1]*Anorm/Bnorm

  return Anorm,Bnorm
end 
#---------------------------------------------------------------------- 

function renormalize_system_touber!(pars)

  ndec    = 5

  f(x,y)  = FXY(x,y,pars.fc0,pars.fcx,pars.fcy)

  f00(y)    = f(0.0,y)
  roots_f   = find_zeros(f00,-10.0,10.0)
  sort!(roots_f)
  α1_2      = roots_f[2]      # Make second root == 1.0

  # A'      = A/Anorm
  Anorm     = α1_2
  

  β1_1      = abs(pars.fcx[1])

  # B       = B/Bnorm
  Bnorm     = Anorm/β1_1 

  # f-parameters

  pars.fc0        = pars.fc0/Anorm
  pars.fcx[1]     = pars.fcx[1]*Bnorm/Anorm

  pars.fcy[1]     = pars.fcy[1]
  pars.fcy[2]     = pars.fcy[2]*Anorm
  pars.fcy[3]     = pars.fcy[3]*(Anorm)^2

  pars.gc0        = pars.gc0/Bnorm
  pars.gcx[1]     = pars.gcx[1]
  pars.gcy[1]     = pars.gcy[1]*Anorm/Bnorm

  return Anorm,Bnorm
end 
#---------------------------------------------------------------------- 
function renormalize_λsystem!(parsS)

  ndec      = 5

  h(x,y)    = FXY(x,y,parsS.fc0,parsS.fcx,parsS.fcy)

  h00(y)    = h(0.0,y)
  roots_f   = find_zeros(h00,-10.0,10.0)
  sort!(roots_f)
  γ2        = roots_f[3]      # Make third root == 1.0

  # λ'      = λ'/λnorm
  λnorm     = γ2

  # h-parameters
  parsS.fc0        = parsS.fc0/λnorm
  parsS.fcx[1]     = parsS.fcx[1]/λnorm

  parsS.fcy[1]     = parsS.fcy[1]
  parsS.fcy[2]     = parsS.fcy[2]*λnorm
  parsS.fcy[3]     = parsS.fcy[3]*(λnorm)^2

  return λnorm
end 
#---------------------------------------------------------------------- 










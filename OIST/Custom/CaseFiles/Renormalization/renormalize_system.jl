function renormalize_system!(pars)

  ndec    = 5

  f(x,y)  = FXY(x,y,pars.fc0,pars.fcx,pars.fcy)

  f00(y)    = f(0.0,y)
  roots_f   = find_zeros(f00,-10.0,10.0)
  sort!(roots_f)
  α1_2      = roots_f[2]      # Make second root == 1.0

  # f-parameters

  pars.fc0        = pars.fc0/α1_2
  pars.fcx[1]     = pars.fcx[1]

  pars.fcy[1]     = pars.fcy[1]
  pars.fcy[2]     = pars.fcy[2]*α1_2
  pars.fcy[3]     = pars.fcy[3]*(α1_2)^2

  pars.gc0        = pars.gc0/α1_2

  return α1_2
end 





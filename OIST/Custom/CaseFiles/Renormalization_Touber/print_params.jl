function print_params(f,g,h,pars,parsS)

  ndec = 5

  # f-parameters
  α1        = pars.fcy[3]
  α3        = pars.fcx[1]

  f00(y)    = f(0.0,y,0.0)
  roots_f   = find_zeros(f00,-10.0,10.0)

  j = 0
  α2 = -999.9
  for i in 1:length(roots_f)
    r = roots_f[i]
    if (r > 1.1)
      α2    = r
    end
  end

  # g-parameters
  β0        = pars.gc0
  β1        = pars.gcy[1]
  β2        = pars.gcx[1]
  R         = sqrt(β1^2 + β2^2)
  ϕ         = atan(β2,β1)

  # h-params
  γ0        = parsS.fc0
  γ1        = parsS.fcx[1]
  
  γ2        = parsS.fcy[3]
  abar      = 0.00
  h0(y)     = h(abar,y)
  hroots_0  = find_zeros(h0,-10.0,10.0)

  @printf "f-params\n"
  @printf "α1 = %.*f\n" ndec α1
  @printf "α2 = %.*f\n" ndec α2
  @printf "α3 = %.*f\n" ndec α3

  @printf "g-params\n"
  @printf "β0   = %.*f\n" ndec β0
  @printf "β1   = %.*f\n" ndec β1
  @printf "β2   = %.*f\n" ndec β2
  @printf "R    = %.*f\n" ndec R
  @printf "ϕ    = %.*f\n" ndec ϕ

  @printf "h-params\n"
  @printf "γ1   = %.*f\n" ndec γ1
  @printf "γ2   = %.*f\n" ndec γ2

  return nothing
end 



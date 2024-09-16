function print_params(f,g,h,pars,parsS)

  ndec = 5

  # f-parameters
  α1_1      = pars.fcy[3]
  β1_1      = pars.fcx[1]

  f00(y)    = f(0.0,y,0.0)
  roots_f   = find_zeros(f00,-10.0,10.0)
  #α1_2      = roots_f[1]
  #α1_3      = roots_f[2]

  nonzeroroots = zeros(Float64,2)
  j = 0
  for i in 1:length(roots_f)
    r = roots_f[i]
    if (abs(r) > 1.0e-6) 
      j = j+1
      nonzeroroots[j] = r
    end
  end
  α1_2      = nonzeroroots[1]
  α1_3      = nonzeroroots[2]

  # g-parameters
  α2_0      = pars.gc0
  α2_1      = pars.gcy[1]
  β2_1      = pars.gcx[1]
  R         = sqrt(α2_1^2 + β2_1^2)
  ϕ         = atan(pars.gcx[1],pars.gcy[1])

  # h-params
  α3_0      = parsS.fc0
  α3_1      = parsS.fcx[1]
  
  γ3_1      = parsS.fcy[3]
  abar      = 0.00
  h0(y)     = h(abar,y)
  hroots_0  = find_zeros(h0,0.001,10.0)
  γ3_2      = hroots_0[1]

  @printf "f-params\n"
  @printf "α1_1 = %.*f\n" ndec α1_1
  @printf "α1_2 = %.*f\n" ndec α1_2
  @printf "α1_3 = %.*f\n" ndec α1_3
  @printf "β1_1 = %.*f\n" ndec β1_1

  @printf "g-params\n"
  @printf "α2_0 = %.*f\n" ndec α2_0
  @printf "α2_1 = %.*f\n" ndec α2_1
  @printf "β2_1 = %.*f\n" ndec β2_1
  @printf "R    = %.*f\n" ndec R
  @printf "ϕ    = %.*f\n" ndec ϕ

  @printf "h-params\n"
  @printf "α3_0 = %.*f\n" ndec α3_0
  @printf "α3_1 = %.*f\n" ndec α3_1
  @printf "γ3_1 = %.*f\n" ndec γ3_1
  @printf "γ3_2 = %.*f\n" ndec γ3_2

  return nothing
end 





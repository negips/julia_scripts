function RenormalizeSystem!(pars)
  # Assuming pars is a mutable structure

  A0              = sqrt(abs(pars.fcy[3]))
  B0              = A0

  @printf "A0: %.5f\n" A0
  @printf "B0: %.5f\n" B0

  # f-system
  pars.fc0        = pars.fc0*A0
  pars.fcy[1]     = pars.fcy[1]
  pars.fcy[2]     = pars.fcy[2]/A0
  pars.fcy[3]     = pars.fcy[3]/(A0^2)

  pars.fcx[1]     = pars.fcx[1]*A0/B0

  # g-system
  pars.gc0        = pars.gc0*B0
  pars.gcy[1]     = pars.gcy[1]*B0/A0
  pars.gcx[1]     = pars.gcx[1]

  return A0,B0
end
#----------------------------------------------------------------------
function RenormalizeλSystem!(pars,A0,B0)

  # h-system
  pars.fcx[1]     = pars.fcx[1]/A0 

  return nothing
end

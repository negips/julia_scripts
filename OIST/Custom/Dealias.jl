function Dealias!(Bcdu,c,du,Geom)

  nn = size(c,2)
  (nd,n1,nel) = size(Geom.gradxd)

  TO = typeof(c[1])

  tmpdu = zeros(TO,n1)
  tmpc  = zeros(TO,n1)

  if nn == 1 # Single column
    for e in 1:nel
      i1 = (e-1)*n1 + 1
      i2 = e*n1

      tmpdx        = Geom.gradxd[:,:,e]*du[i1:i2]       # Du
      tmpc         = Geom.intpm1d*c[i1:i2]              # Int*v
      tmpc         = Geom.bm1d[:,e].*tmpc.*tmpdx        # B*(Int*v)*Du

      Bcdu[i1:i2]  = Geom.intpm1d'*tmpc
    end   

  else
   
    for e in 1:nel
      tmpdx        = Geom.gradxd[:,:,e]*du[:,e]         # Du
      tmpc         = Geom.intpm1d*c[:,e]                # Int*v
      tmpc         = Geom.bm1d[:,e].*tmpc.*tmpdx        # B*(Int*v)*Du

      Bcdu[:,e]    = Geom.intpm1d'*tmpc
    end

  end  

end  

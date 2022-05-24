#     Julia Port for speclib.f (Original Author: Einar Malvin Ronquist)
#     From Nek5000.
#     Author:     Prabal Negi
#

      module SpecLib

      export jacobf

#      function __init__()
#
#        if rank == 0
#          println("Initialied MPI in Module JNek_IO")
#        end  
#
#        return nothing
#      end 
#----------------------------------------------------------------------

      function jacobf(n::Int,α::Float64,β::Float64,x::Float64)

#     Computes the Jacobi polynomial (p) and its derivative (pd)
#     of degree N at X.

      p::Float64 = 0.0
      pd::Float64 = 0.0

      if (n==0)
        return p,pd
      end  

      apb  = α+β
      p  = 1.
      pd = 0.

      polyl = p
      pderl = pd
      p     = (α-β+(apb+2.)*x)/2.0
      pd    = (apb+2.)/2.0

      for k in 2:n
         dk     = k+0.0
         a1     = 2.0*dk*(dk+apb)*(2.0*dk+apb-2.0)
         a2     = (2.0*dk+apb-1.0)*(α^2.0-β^2.0)
         b3     = (2.0*dk+apb-2.0)
         a3     = b3*(b3+1.0)*(b3+2.0)
         a4     = 2.0*(dk+α-1.0)*(dk+β-1.0)*(2.0*dk+apb)
         polyn  = ((a2+a3*x)*p-a4*polyl)/a1
         pdern  = ((a2+a3*x)*pd-a4*pderl+a3*p)/a1
         psave  = polyl
         pdsave = pderl
         polyl  = p
         p      = polyn
         pderl  = pd
         pd     = pdern
      end   
#      polym1 = polyl
#      pderm1 = pderl
#      polym2 = psave
#      pderm2 = pdsave

      return p,pd
      end
#---------------------------------------------------------------------- 



      end   # Module

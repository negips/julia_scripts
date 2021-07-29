      include("BulgeChase.jl")
      Hes = copy(Hold)
      VV  = Vold
      kk = LKryl

      OH = Hes[1:kk,1:kk]
      Hc = copy(OH)
      μ,nμ  = ArnGetLowerShifts(Hc,EKryl)

      Hf = copy(Hc)
      for j in 1:1
        global Hf,Hp,Qf1
        local Qp1,nμ
        nμ    = 1
        Hp,Qp1 = CreateBulge(Hf,μ[j:j],nμ);
        Hf,Qf1 = ChaseBulgeDown(Hp,nμ);
      end  

      k  = 1
      k1 = 2
      k2 = 3
      x = Hp[:,k]
      qi,w,τ = CreateReflectorZerosSub(x,2,3,kk)
#      Hsub     = Hp[k:k2,k:k2]
#      a = qi*Hsub
#      b = a*qi'
#      Hp[k1:k2,k1:k2] = 1.0*Hsub

      println("done")




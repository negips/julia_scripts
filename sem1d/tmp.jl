      include("BulgeChase.jl")

      Hes = copy(Hold)
      VV  = Vold
      kk = LKryl

      OH = Hes[1:kk,1:kk]
      Hc = copy(OH)
      μ,nμ  = ArnGetLowerShifts(Hc,EKryl)

      A   = deepcopy(OH)
      B   = deepcopy(OH)

      rw,cl = size(OH)

      Q   = Matrix{typeof(OH[1,1])}(1.0I,rw,cl)
      T   = Matrix{typeof(OH[1,1])}(1.0I,rw,cl)      # tmp
      Qi  = Matrix{typeof(OH[1,1])}(1.0I,rw,cl)

      tol = 1.0e-12
      Hf = copy(Hc)
      tμ = nμ
      for j in 1:tμ
        global Hf,Hp,Qf1,Qp1,Qi,Hfold,τ
        local nλ
        global wi
        local A
        nλ    = 1
        Hp,Qp1 = CreateBulge(Hf,μ[j:j],nλ);
        Hf,Qf1 = ChaseBulgeDown(Hp,nλ);

      end  

      qq = Qf1*Qp1;

      tt = LKryl-tμ
      βk = Hf[tt+1,tt]
      println("βk=$βk; done")




      include("BulgeChase.jl")

      Hes = copy(Hold)
      VV  = Vold
      kk  = LKryl

      EKryl = 10 #37

      OH = Hes[1:kk,1:kk]
      Hc = copy(OH)
      μ,nμ  = ArnGetLowerShifts(Hc,EKryl)

      A   = deepcopy(OH)
      B   = deepcopy(OH)

      rw,cl = size(OH)

      Q   = Matrix{eltype(OH)}(1.0I,rw,cl)
      T   = Matrix{eltype(OH)}(1.0I,rw,cl)      # tmp
      Qi  = Matrix{eltype(OH)}(1.0I,rw,cl)

      tol = 1.0e-12
      Hf = copy(OH)
      tμ = EKryl
      for j in 1:tμ
        for k in 1:5
          global Hf,Hp,Qf1,Qp1,Qi,Hfold,τ
          local nλ
          global wi
          local A, λ
          nλ    = 1
          λ     = μ[j:j]
          Hp,Qp1 = CreateBulge(Hf,λ,nλ);
          Hf,Qf1 = ChaseBulgeDown1(Hp,λ,nλ);
        end               
      end  

      qq = Qf1*Qp1;

      tt = LKryl-tμ
      βk1 = Hf[tt+1,tt]


      nλ = EKryl
      λ  = μ[1:nλ]
      Hfr,Qfr = FrancisSeq(Hc,λ,nλ)
      tt2 = LKryl-nλ
      βk2 = Hfr[tt2+1,tt2] 

      nλ = EKryl
      λ  = μ[1:nλ]
      Hqr,Qqr = ExplicitShiftedQR(Hc,λ,nλ,5)
      tt3 = LKryl-nλ
      βk3 = Hqr[tt2+1,tt2] 

      nλ = EKryl
      λ  = μ[1:nλ]
      Hqr2,Qqr2 = ExplicitShiftedQR(Hqr,λ,nλ,ngs)
      tt4 = LKryl-nλ
      βk4 = Hqr2[tt2+1,tt2] 


      println("βk1 (manual)    =$βk1; done")
      println("βk2 (Francis)   =$βk2; done")
      println("βk3 (Explicit)  =$βk3; done")
      println("βk4 (Explicit2) =$βk4; done")




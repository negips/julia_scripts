function AssembleMatrixLesshafft(c0,cnv,wlp,xm1,bm1,Basis,lx1,nel)

#     Building the Complex Ginzburg Landau model problem from
#     Lutz Lesshafft (2018) Artificial eigenmodes in truncated flow domains
#     
#     ∂ψ/∂t     = -U∂ψ/∂x + μ(x)ψ + γ∂∂ψ/∂x∂x + f(x)ψ(s)
#     f(x)ψ(xₛ) = c0*exp(-((x-xₐ)/0.1)²)*ψ(s)
#     xₛ        = 39.0
#     a         = 1.0
#     c0        = Parametrically varied.
#
#     μ(x)      = U²/8*(1 - x/20)
#     U         = 6
#
#     γ         = 1. - i

#      VT  = Float64
      VT  = Complex

      II  = Matrix(1.0I,lx1,lx1)
#      Iv  = Vector(1.0,lx1)
      U   = 6.0
      γ   = 1.0 - 1.0*im

      dof = nel*(lx1-1) + 1;

      A    = zeros(VT,dof,dof);
      B    = zeros(Float64,dof,dof);
      Binv = zeros(Float64,lx1,nel);

      OP   = zeros(VT,lx1,lx1,nel);

      xa  = 1.0         # Feedback destination point
      xs  = 39.0        # Feedback source point
      b   = 0.1         # Exponential drop off rate for feedback

#     Find the element containing xₛ
      ar   = argmin((xm1 .- xs).^2)
#      r    = ar.I[1]
      e    = ar.I[2]          # Element which has the point xₛ

#     Map xₛ to reference coordinates ξ in [-1,1]
      ξ    = map_to_canonical(xs,xm1[1,e],xm1[lx1,e],Basis)
      Ixs  = zeros(Float64,1,lx1)
      interpolation_matrix!(Ixs,ξ,zgm1,Basis.baryweights)

#     Location of Ixs in the global matrix (Columns)
      je1 = (e-1)*(lx1-1) + 1;
      je2 = je1 + lx1 - 1

      for i in 1:nel
        j1 = (i-1)*(lx1-1) + 1;
        j2 = j1 + lx1 -1;

#       Mass matrix
        B[j1:j2,j1:j2] = B[j1:j2,j1:j2] + diagm(bm1[:,i])

     
#       Standard feedback term
        μ   = (U*U/8.0)*(1.0 .- xm1[:,i]./20.0)
        bμ  = bm1[:,i].*μ
        Mμ  = diagm(bμ)
#        display(Mμ)

#       Sub matrix 
        subm = -U.*cnv[:,:,i] + Mμ + γ.*wlp[:,:,i];
        OP[:,:,i] = subm
        A[j1:j2,j1:j2] = A[j1:j2,j1:j2] + subm;

#       Long range Feedback Term
        expx = c0.*exp.(-((xm1[:,i] .- xa)/b).^2)
        Fmax = maximum(expx)
        if Fmax > 1.0e-14
          bmex = bm1[:,i].*expx
          Feed = bmex*Ixs
          A[j1:j2,je1:je2] = A[j1:j2,je1:je2] + Feed
          OP[:,:,i] = OP[:,:,i] + Feed
        end  
      end

#     Binv Vector (Since its a Diagonal matrix)
      for i in 1:nel
        j1 = (i-1)*(lx1-1) + 1;
        j2 = j1 + lx1 -1;
        b  = diag(B[j1:j2,j1:j2])
        Binv[:,i] = 1.0./b
      end 

      println("Matrix Built")
      return A, B, OP, Binv
end








"""
    AssembleMatrixCRD(U,γ,cnv,wlp,xm1,bm1,Basis,lx1,nel,prec)

     Building the Matrices of a Convection Reaction Diffusion model:
     
     ∂ψ/∂t     = -U∂ψ/∂x + μ(x)ψ + γ∂∂ψ/∂x∂x

     μ(x)      = 0*(x)

     γ     - Diffusion Coefficient
     Conv  - Convection term:                  -Udψ/dx
     Src   - Spatially varying source term:    μ(x)ψ
     Lap   - Weak Laplacian:                   γd²ψ/dx²
     B     - Mass matrix:
     A     - Combined operator:    -Udψ/dx + μ(x)ψ + γd²ψ/dx² + F(x,s)ψ

"""
function AssembleMatrixCRD(U,γ,cnv,wlp,xm1,bm1,Basis,lx1,nel,prec,ifsparse)

      VT  = prec

      II  = Matrix{VT}(I,lx1,lx1)

      zro   = prec(0)
      one   = prec(1)
      half  = prec(0.5)

      dof = nel*lx1;

      if ifsparse
        A    = spzeros(VT,dof,dof);
        B    = zeros(VT,dof);
        Conv = spzeros(VT,dof,dof)
        Src  = spzeros(VT,dof,dof)
        Lap  = spzeros(VT,dof,dof)
      else
        A    = zeros(VT,dof,dof);
        B    = zeros(VT,dof);
        Conv = zeros(VT,dof,dof)
        Src  = zeros(VT,dof,dof)
        Lap  = zeros(VT,dof,dof)
      end  

      OP   = zeros(VT,lx1,lx1,nel);

#     Mass matrix
      B = bm1[:]

      for i in 1:nel
        j1 = (i-1)*lx1 + 1;
        j2 = i*lx1;

        μ   = ones(VT,length(lx1))
        bμ  = bm1[:,i].*μ
        Mμ  = diagm(bμ)

        Conv[j1:j2,j1:j2] = -U.*cnv[:,:,i]
        Src[j1:j2,j1:j2]  = Mμ
        Lap[j1:j2,j1:j2]  = γ.*wlp[:,:,i]

#       Sub matrix 
        subm = -U.*cnv[:,:,i] + Mμ + γ.*wlp[:,:,i];
        OP[:,:,i] = subm
        A[j1:j2,j1:j2] = A[j1:j2,j1:j2] + subm;
      end

      println("Direct Sparse Matrices Built")
      return A, B, OP, Conv, Src, Lap
end
#---------------------------------------------------------------------- 

function AssembleFilter(fil,M2N,N2M,lx1,nel,prec,ifsparse)

#     M2N - Modal to Local Transformation
#     N2M - Nodal to Modal Transformation

      VT  = prec

      II  = Matrix{VT}(I,lx1,lx1)

      dof = nel*lx1;

      if ifsparse
        Fil  = spzeros(VT,dof,dof)
      else
        Fil  = zeros(VT,dof,dof);
      end  

      OP  = zeros(VT,lx1,lx1,nel)

      for i in 1:nel
        j1 = (i-1)*lx1 + 1;
        j2 = i*lx1;

        subm              = M2N*diagm(fil)*N2M
        Fil[j1:j2,j1:j2]  = subm 
        OP[:,:,i]         = subm
      end

      println("Filter Matrix Built")
      return Fil, OP
end
#---------------------------------------------------------------------- 

function BuildFilter(M2N,N2M,lx1,nel,prec,ifsparse)

      VT  = prec

      f   = zeros(VT,lx1)
      χ   = 1.0

      f[lx1]   = 1.0
#      f[lx1-1] = 0.5
#      f[lx1-2] = 0.25

      f        = χ*f

      Fil,OPf = AssembleFilter(f,M2N,N2M,lx1,nel,prec,ifsparse)

      return Fil,OPf
end
#---------------------------------------------------------------------- 
function AssembleInterpolation(intp,nel,ifsparse)

#     intp - Local interpolation Matrix

      VT  = eltype(intp[1])

      lxd,lx1 = size(intp)

      dof1 = nel*lx1;
      dof2 = nel*lxd;

      if ifsparse
        Intp  = spzeros(VT,dof2,dof1)
      else
        Intp  = zeros(VT,dof2,dof1);
      end  

      for i in 1:nel
        j1 = (i-1)*lxd + 1;
        j2 = i*lxd;

        j3 = (i-1)*lx1 + 1;
        j4 = i*lx1;

        Intp[j1:j2,j3:j4]  = intp 
      end

      println("Interpolation Matrix Built")
      return Intp
end
#---------------------------------------------------------------------- 







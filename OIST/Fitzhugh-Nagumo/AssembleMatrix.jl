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
    
##       Standard source term
#        if (whichsrc==1)
#          if (i==1)
#            @printf "Using Linear Decrease x/%3.2f for source term.\n" cx0 
#          end  
#          μ   = xm1[:,i]
#        elseif (whichsrc==2) 
#          if i==1
#            println("Using 0.5*tanh(x/20) for source term.")
#          end  
#          μ   = half*tanh.(xm1[:,i]./twnty)
#        else
#          if i==1
#            println("whichsrc=$whichsrc not defined")
#          end  
#        end  

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
# function AssembleAdjointLesshafftSparse(U,γ,c0,cx0,whichsrc,cnv,wlp,xm1,bm1,Basis,lx1,nel,prec)
# 
# #     Building the Complex Ginzburg Landau model problem from
# #     Lutz Lesshafft (2018) Artificial eigenmodes in truncated flow domains
# #     
# #     ∂ψ/∂t     = -U∂ψ/∂x + μ(x)ψ + γ∂∂ψ/∂x∂x + f(x)ψ(s)
# #     f(x)ψ(xs) = c0*exp(-((x-xₐ)/0.1)²)*ψ(s)
# #     xs        = 39.0
# #     xa        = 1.0
# #     c0        = Parametrically varied.
# #
# #     cx0       = Point where system is no longer convectively unstable
# #     μ(x)      ==      1:  U²/8*(1 - x/cx0)
# #                       2:  U²/8*(1 - 0.5*tanh(x))
# #     U         = 6
# #
# #     γ         = 1. - i
# 
# #     Conv  - Convection term:                  -Udψ/dx
# #     Src   - Spatially varying source term:    μ(x)ψ
# #     Lap   - Weak Laplacian:                   γd²ψ/dx²
# #     B     - Mass matrix:
# #     Fd    - Feedback Matrix:                  F(x,s)ψ
# #     A     - Combined operator:    -Udψ/dx + μ(x)ψ + γd²ψ/dx² + F(x,s)ψ
# 
#       VT  = Complex{prec}
# 
#       II  = Matrix{VT}(I,lx1,lx1)
# 
#       xa    = prec(1)           # Feedback destination point
#       xs    = prec(39)          # Feedback source point
#       b     = prec(0.1)         # Exponential drop off rate for feedback
# 
#       fact  = prec(100)
#       zro   = prec(0)
#       one   = prec(1)
#       eight = prec(8)
#       twnty = prec(20)
#       half  = prec(0.5)
# 
#       cutoff = fact*eps(prec)
# 
#       if (abs(c0)>cutoff)
#         @printf "Feedback not implemented for Adjoint Equations yet.\n"
#         @printf "Setting c0=0.0\n"
#         c0 = zro
#       end  
# 
#       dof = nel*lx1;
# 
#       A    = spzeros(VT,dof,dof);
#       B    = zeros(VT,dof);
#       Conv = spzeros(VT,dof,dof)
#       Src  = spzeros(VT,dof,dof)
#       Lap  = spzeros(VT,dof,dof)
#       Fd   = spzeros(VT,dof,dof)
#       Binv = zeros(VT,lx1,nel);
# 
#       OP   = zeros(VT,lx1,lx1,nel);
# 
# 
# #     Find the element containing xs
#       ar   = argmin((xm1 .- xs).^2)
#       e    = ar.I[2]          # Element which has the point xs
# 
# #     Map xs to reference coordinates ξ in [-1,1]
#       ξ    = map_to_canonical(xs,xm1[1,e],xm1[lx1,e],Basis)
#       Ixs  = zeros(VT,1,lx1)
#       interpolation_matrix!(Ixs,ξ,zgm1,Basis.baryweights)
# 
# #     Location of Ixs in the global matrix (Columns)
#       je1 = (e-1)*lx1 + 1;
#       je2 = e*lx1
# 
# #     Mass matrix
#       B     = bm1[:]
#       γ2    = γ'
# 
#       for i in 1:nel
#         j1 = (i-1)*lx1 + 1;
#         j2 = i*lx1;
#     
# #       Standard source term
#         if (whichsrc==1)
#           if i==1
#             @printf "Using Linear Decrease x/%3.2f for source term.\n" cx0 
#           end  
#           μ   = (U*U/eight)*(one .- xm1[:,i]./cx0)
#         elseif (whichsrc==2) 
#           if i==1
#             println("Using 0.5*tanh(x/20) for source term.")
#           end  
#           μ   = (U*U/eight)*(one .- half*tanh.(xm1[:,i]./twnty))
#         else
#           if i==1
#             println("whichsrc=$whichsrc not defined")
#           end  
#         end  
#        
#         bμ  = bm1[:,i].*μ
#         Mμ  = diagm(bμ)
# 
# #       Sign for convection term changes in Adjoint        
#         Conv[j1:j2,j1:j2] = U.*cnv[:,:,i]
#         Src[j1:j2,j1:j2]  = Mμ
#         Lap[j1:j2,j1:j2]  = γ2.*wlp[:,:,i]
#         if i==nel
#           Lap[j2,j2]      = Lap[j2,j2] - U/γ2
#         end 
# 
# #       Sub matrix 
#         subm = U.*cnv[:,:,i] + Mμ + γ2.*wlp[:,:,i];
#         OP[:,:,i] = subm
#         A[j1:j2,j1:j2] = A[j1:j2,j1:j2] + subm;
# 
# #       Long range Feedback Term
#         expx = c0.*exp.(-((xm1[:,i] .- xa)/b).^2)
#         Fmax = maximum(abs.(expx))
#         if Fmax > cutoff
#           bmex = bm1[:,i].*expx
#           Feed = c0.*bmex*Ixs
#           A[j1:j2,je1:je2] = A[j1:j2,je1:je2] + Feed
# 
#           Fd[j1:j2,je1:je2] = copy(Feed)
#         end  
#       end
# 
#       println("Adjoint Sparse Matrices Built")
#       return A, B, OP, Conv, Src, Lap, Fd
# end
# #---------------------------------------------------------------------- 

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

      println("Filter Matricx Built")
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






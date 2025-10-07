function GinzburgLandauLocalOP(δ::Vector{T},Inp::SEMInput,GeoM::SEMGeoMat{GT},B::NodalBasis) where {T<:Number,GT<:Number}

#     Building the Complex Ginzburg Landau model problem from

#     U           = δ[1]
#     μ0          = δ[2]
#     μ           = δ[3]
#     γ           = δ[4]
     
#     ∂ψ/∂t       = δ[1]∂ψ/∂x + δ[2]ψ + δ[3]xψ + δ[4]∂(∂ψ/∂x)/∂x

#     Conv  - Convection term:                  δ[1]dψ/dx
#     Src   - Spatially varying source term:    (δ[2] + δ[3]x)ψ
#     Lap   - Weak Laplacian:                   δ[4]d²ψ/dx²
#     B     - Mass matrix:
#     A     - Combined Linear operator:    Lψ = (δ[1]d/dx + δ[2]I + δ[3]x + δ[4]d²/dx²)ψ

      Prec  = Complex{Inp.Dtype}
      lx1   = Inp.lx1
      nel   = Inp.nel

      II    = Matrix{Prec}(I,lx1,lx1)

      zro   = Prec(0)
      one   = Prec(1)
      eight = Prec(8)
      twnty = Prec(20)
      half  = Prec(0.5)
      hund  = Prec(100)

      OP   = zeros(Prec,lx1,lx1,nel);

      for i in 1:nel
        j1 = (i-1)*lx1 + 1;
        j2 = i*lx1;
    
        # Standard source term
        μ         = (δ[2] .+ GeoM.xm1[:,i]*δ[3])
        subm      = δ[1].*GeoM.Convd[:,:,i] .+ diagm(GeoM.bm1[:,i].*μ) .+ δ[4]*GeoM.WLap[:,:,i]
        OP[:,:,i] = subm

      end

      println("Direct Local Matrices Built")
      return OP
end
#---------------------------------------------------------------------- 

function GinzburgLandauSparse(δ::Vector{T},Inp::SEMInput,GeoM::SEMGeoMat{GT},B::NodalBasis) where {T<:Number,GT<:Number}

#     Building the Complex Ginzburg Landau model problem

#     U           =-δ[1]
#     μ0          = δ[2]
#     μ           = δ[3]
#     γ           = δ[4]
     
#     ∂ψ/∂t       = δ[1]∂ψ/∂x + δ[2]ψ + δ[3]xψ + δ[4]∂(∂ψ/∂x)/∂x

#     Conv  - Convection term:                  δ[1]dψ/dx
#     Src   - Spatially varying source term:    (δ[2] + δ[3]x)ψ
#     Lap   - Weak Laplacian:                   δ[4]d²ψ/dx²
#     B     - Mass matrix:
#     A     - Combined Linear operator:    Lψ = (δ[1]d/dx + δ[2]I + δ[3]x + δ[4]d²/dx²)ψ

      Prec  = Complex{Inp.Dtype}
      lx1   = Inp.lx1
      nel   = Inp.nel

      II    = Matrix{Prec}(I,lx1,lx1)

      zro   = Prec(0)
      one   = Prec(1)
      eight = Prec(8)
      twnty = Prec(20)
      half  = Prec(0.5)
      hund  = Prec(100)


      dof  = nel*lx1;

      A    = spzeros(Prec,dof,dof);
      Conv = spzeros(Prec,dof,dof)
      Src  = spzeros(Prec,dof,dof)
      Lap  = spzeros(Prec,dof,dof)

      OP   = zeros(Prec,lx1,lx1,nel);

      # Mass matrix
      B    = GeoM.bm1[:]

      for i in 1:nel
        j1 = (i-1)*lx1 + 1;
        j2 = i*lx1;
    
        # Standard source term
        μ   = (δ[2] .+ GeoM.xm1[:,i]*δ[3])

        bμ  = GeoM.bm1[:,i].*μ
        Mμ  = diagm(bμ)

        Conv[j1:j2,j1:j2] = δ[1]*GeoM.Convd[:,:,i]
        Src[j1:j2,j1:j2]  = diagm(GeoM.bm1[:,i].*μ)
        Lap[j1:j2,j1:j2]  = δ[4]*GeoM.WLap[:,:,i]

        # Sub matrix 
        subm            = δ[1].*GeoM.Convd[:,:,i] + diagm(GeoM.bm1[:,i].*μ) + δ[4]*GeoM.WLap[:,:,i]
        OP[:,:,i]       = subm
        A[j1:j2,j1:j2]  = A[j1:j2,j1:j2] + subm

      end

      println("Direct Sparse Matrices Built")
      return A, B, OP, Conv, Src, Lap
end
#---------------------------------------------------------------------- 
function AssembleAdjointGLSparse(U,γ,μ0,μx,whichsrc,gradx,cnv,wlp,xm1,bm1,Basis,lx1,nel,prec)

#     Building the Complex Ginzburg Landau model problem from
#     Lutz Lesshafft (2018) Artificial eigenmodes in truncated flow domains
#     
#     ∂ψ/∂t     = -U∂ψ/∂x + μ(x)ψ + γ∂∂ψ/∂x∂x + f(x)ψ(s)
#     f(x)ψ(xs) = c0*exp(-((x-xₐ)/0.1)²)*ψ(s)
#     xs        = 39.0
#     xa        = 1.0
#     c0        = Parametrically varied.
#
#     cx0       = Point where system is no longer convectively unstable
#     μ(x)      ==      1:  U²/8*(1 - x/cx0)
#                       2:  U²/8*(1 - 0.5*tanh(x))
#     U         = 6
#
#     γ         = 1. - i

#     Conv  - Convection term:                  -Udψ/dx
#     Src   - Spatially varying source term:    μ(x)ψ
#     Lap   - Weak Laplacian:                   γd²ψ/dx²
#     B     - Mass matrix:
#     Fd    - Feedback Matrix:                  F(x,s)ψ
#     A     - Combined operator:    -Udψ/dx + μ(x)ψ + γd²ψ/dx² + F(x,s)ψ

      VT  = Complex{prec}

      II  = Matrix{VT}(I,lx1,lx1)

      xa    = prec(1)           # Feedback destination point
      xs    = prec(39)          # Feedback source point
      b     = prec(0.1)         # Exponential drop off rate for feedback
      c0    = prec(0.0)

      fact  = prec(100)
      zro   = prec(0)
      one   = prec(1)
      eight = prec(8)
      twnty = prec(20)
      half  = prec(0.5)

      cutoff = fact*eps(prec)

      # if (abs(c0)>cutoff)
      #   @printf "Feedback not implemented for Adjoint Equations yet.\n"
      #   @printf "Setting c0=0.0\n"
      #   c0 = zro
      # end  

      dof = nel*lx1;

      A    = spzeros(VT,dof,dof);
      B    = zeros(VT,dof);
      Conv = spzeros(VT,dof,dof)
      Src  = spzeros(VT,dof,dof)
      Lap  = spzeros(VT,dof,dof)
      SLap = spzeros(VT,dof,dof)
      Fd   = spzeros(VT,dof,dof)
      Binv = zeros(VT,lx1,nel);

      OP   = zeros(VT,lx1,lx1,nel);


#     Find the element containing xs
      ar   = argmin((xm1 .- xs).^2)
      e    = ar.I[2]          # Element which has the point xs

#     Map xs to reference coordinates ξ in [-1,1]
      ξ    = map_to_canonical(xs,xm1[1,e],xm1[lx1,e],Basis)
      Ixs  = zeros(VT,1,lx1)
      interpolation_matrix!(Ixs,ξ,zgm1,Basis.baryweights)

#     Location of Ixs in the global matrix (Columns)
      je1 = (e-1)*lx1 + 1;
      je2 = e*lx1

#     Mass matrix
      B     = bm1[:]
      γ2    = γ'

      for i in 1:nel
        j1 = (i-1)*lx1 + 1;
        j2 = i*lx1;
    
#       Standard source term
        if (whichsrc==1)
          # if i==1
          #   @printf "Using Linear Decrease x/%4.3f for source term.\n" abs(μx)
          # end  
          μ   = (μ0' .+ xm1[:,i]*μx')
        else
          if i==1
            println("whichsrc=$whichsrc not defined")
          end  
        end  
       
        bμ  = bm1[:,i].*μ
        Mμ  = diagm(bμ)

        # Sign for convection term changes in Adjoint        
        Conv[j1:j2,j1:j2] = U.*cnv[:,:,i]
        Src[j1:j2,j1:j2]  = Mμ
        Lap[j1:j2,j1:j2]  = γ2.*wlp[:,:,i]
        if i==nel
          Lap[j2,j2]      = Lap[j2,j2] - U/γ2
        end 
        SLap[j1:j2,j1:j2] = gradx[:,:,i]*gradx[:,:,i]

        # Sub matrix 
        subm = U.*cnv[:,:,i] + Mμ + γ2.*wlp[:,:,i];
        OP[:,:,i] = subm
        A[j1:j2,j1:j2] = A[j1:j2,j1:j2] + subm;

        # Long range Feedback Term
        expx = c0.*exp.(-((xm1[:,i] .- xa)/b).^2)
        Fmax = maximum(abs.(expx))
        if Fmax > cutoff
          bmex = bm1[:,i].*expx
          Feed = c0.*bmex*Ixs
          A[j1:j2,je1:je2] = A[j1:j2,je1:je2] + Feed

          Fd[j1:j2,je1:je2] = copy(Feed)
        end  
      end

      println("Adjoint Sparse Matrices Built")
      return A, B, OP, Conv, Src, Lap, Fd, SLap
end
#---------------------------------------------------------------------- 
function AdjointGinzburgLandauSparse(δ::Vector{T},Inp::SEMInput,GeoM::SEMGeoMat{GT},B::NodalBasis) where {T<:Number,GT<:Number}


#     Building the Adjoint Complex Ginzburg Landau model problem

#     U           =-δ[1]
#     μ0          = δ[2]
#     μ           = δ[3]
#     γ           = δ[4]
     
#     ∂ψ*/∂t       = -δ[1]'∂ψ*/∂x + δ[2]'ψ* + δ[3]'xψ* + δ[4]'∂(∂ψ*/∂x)/∂x

#     Conv  - Convection term:                  -δ[1]'dψ*/dx
#     Src   - Spatially varying source term:    (δ[2]' + δ[3]'x)ψ*
#     Lap   - Weak Laplacian:                   δ[4]'d²ψ*/dx²
#     B     - Mass matrix:
#     A     - Combined Linear operator:    Lψ = (-δ[1]'d/dx + δ[2]'I + δ[3]'x + δ[4]'d²/dx²)ψ*

      Prec  = Complex{Inp.Dtype}
      lx1   = Inp.lx1
      nel   = Inp.nel

      II    = Matrix{Prec}(I,lx1,lx1)

      zro   = Prec(0)
      one   = Prec(1)
      eight = Prec(8)
      twnty = Prec(20)
      half  = Prec(0.5)
      hund  = Prec(100)


      dof   = nel*lx1;

      A     = spzeros(Prec,dof,dof)
      Conv  = spzeros(Prec,dof,dof)
      Src   = spzeros(Prec,dof,dof)
      Lap   = spzeros(Prec,dof,dof)

      OP    = zeros(Prec,lx1,lx1,nel)

      # Mass matrix
      B     = GeoM.bm1[:]

      for i in 1:nel
        j1  = (i-1)*lx1 + 1;
        j2  = i*lx1;
    
        # Standard source term
        μ   = (δ[2]' .+ (δ[3]')*GeoM.xm1[:,i])

        # Sign for convection term changes in Adjoint        
        Conv[j1:j2,j1:j2] = -(δ[1]').*GeoM.Convd[:,:,i]
        Src[j1:j2,j1:j2]  =  diagm(GeoM.bm1[:,i].*μ)
        Lap[j1:j2,j1:j2]  =  (δ[4]').*GeoM.WLap[:,:,i]

        # Boundary Term (for Adjoint)
        if i==1
          Lap[j1,j1]      = Lap[j1,j1] + (-δ[1]/δ[4])'
        elseif i==nel
          Lap[j2,j2]      = Lap[j2,j2] - (-δ[1]/δ[4])'
        end 

        # Sub matrix 
        subm            = -(δ[1]').*GeoM.Convd[:,:,i] + diagm(GeoM.bm1[:,i].*μ) + (δ[4]').*GeoM.WLap[:,:,i]
        OP[:,:,i]       = subm
        A[j1:j2,j1:j2]  = A[j1:j2,j1:j2] + subm;

      end

      println("Adjoint Sparse Matrices Built")
      return A, B, OP, Conv, Src, Lap
end
#---------------------------------------------------------------------- 
function GLSetBC!(A::AbstractMatrix{T},lbc::Bool,rbc::Bool,ifperiodic::Bool) where {T<:Number}

  r,c = size(A)

  if (ifperiodic)
    # do nothing
  else
    if lbc
      A[1,:] = zeros(T,1,c)
      A[1,1] = T(1)
    end
    if rbc
      A[r,:] = zeros(T,1,c)
      A[r,c] = T(1)
    end
  end  
  return nothing
end
#---------------------------------------------------------------------- 
function SEM_SetBC!(v::AbstractVector{T},lbc::Bool,rbc::Bool) where {T<:Number}


  if lbc
    v[1]   = T(0)
  end
  if rbc
    v[end] = T(0)
  end

  return nothing
end
#---------------------------------------------------------------------- 
function GLAnalyticalSpectra(δ::Vector{T}) where {T<:Number}

  # Analytical Eigenvalues on a semi-infinite domain
  
  # Zeros of the Airy function
  ω1        = find_zero(airyai,(-3.0,-0.0))
  ω2        = find_zero(airyai,(-5.0,-3.0))
  ω3        = find_zero(airyai,(-6.0,-5.0))
  ω4        = find_zero(airyai,(-7.0,-6.0))
  ω5        = find_zero(airyai,(-8.0,-7.0))
  ω6        = find_zero(airyai,(-9.5,-8.0))
  ω7        = find_zero(airyai,(-10.5,-9.5))
  ω8        = find_zero(airyai,(-11.8,-10.5))
  ω9        = find_zero(airyai,(-12.0,-11.8))
  ω10       = find_zero(airyai,(-12.9,-12.0))
  ω11       = find_zero(airyai,(-13.8,-12.9))
  ω12       = find_zero(airyai,(-14.8,-13.8))
  ω13       = find_zero(airyai,(-15.8,-14.8))
  ω14       = find_zero(airyai,(-16.8,-15.8))
  ω15       = find_zero(airyai,(-17.5,-16.8))
  ω         = [ω1, ω2, ω3, ω4, ω5, ω6, ω7, ω8, ω9, ω10, ω11, ω12, ω13, ω14, ω15]

  U         = -δ[1]
  μ0        =  δ[2]
  μx        =  δ[3]
  γ         =  δ[4]

  Ω         = (μ0 .- U*U/(4.0*γ) .+ (γ*μx*μx)^(1.0/3.0)*ω)

  return Ω
end



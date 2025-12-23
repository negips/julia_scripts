# Calculate the Center-Manifold Asymptotic vectors

include("../Module_CenterManifold/CenterManifold.jl")

#---------------------------------------------------------------------- 
function GLStdAsympNLTerm(n0::Int,Ord0::Int,Nc::Int,MV1::AbstractMatrix{T},MV2::AbstractMatrix{T},MV3::AbstractMatrix{T},Ord1::Int,Ord2::Int,Ord3::Int,δ5) where {T}

  ind0 = CenterManifold.GetPolynomialIndices(n0,Ord0,Nc) .+ 1

  Nt1  = CenterManifold.NInteractionTerms(Ord1,Nc)
  Nt2  = CenterManifold.NInteractionTerms(Ord2,Nc)
  Nt3  = CenterManifold.NInteractionTerms(Ord3,Nc)

  N,c  = size(MV1)
  Nby2 = Int(N/2)

  h    = zeros(T,N)

  topi = 1:Nby2
  boti = Nby2+1:N

  if (Ord1 + Ord2 + Ord3 == Ord0) && Ord0>=3

    # The standard nonlinearity δ5*|A*A|A and δ5'|A*A|A* only starts at third-order.
    for i in 1:Nt1
      ind1 = CenterManifold.GetPolynomialIndices(i,Ord1,Nc) .+ 1
      for j in 1:Nt2
        ind2 = CenterManifold.GetPolynomialIndices(j,Ord2,Nc) .+ 1
        for k in 1:Nt3
          ind3 = CenterManifold.GetPolynomialIndices(k,Ord3,Nc) .+ 1

          ind_total = [ind1[:];ind2[:];ind3[:]]

          sort!(ind_total)

          if (ind_total == ind0)    # Matching polynomials
            # println("Index Match: $ind1, $ind2, $ind3: $ind0")
            h[topi] .= h[topi] .+ (δ5 )*MV1[boti,i].*MV2[topi,j].*MV3[topi,k]
            h[boti] .= h[boti] .+ (δ5')*MV1[topi,i].*MV2[boti,j].*MV3[boti,k]
          end            
        end       # k
      end         # j
    end           # i
  end             # if (Ord1 + Ord2 + Ord3 == Ord0) && Ord0>=3

  return h
end
#---------------------------------------------------------------------- 
function GLLapNLTerm(n0::Int,Ord0::Int,Nc::Int,MV1::AbstractMatrix{T},Ord1::Int,mode::Int,modeC::Int,Lap::AbstractMatrix{S},LapC::AbstractMatrix{S}) where {T,S}

  # mode    - Mode No of δ4
  # modeC   - Mode No of δ4' (conjugate)

  ind0 = CenterManifold.GetPolynomialIndices(n0,Ord0,Nc) .+ 1

  Nt1  = CenterManifold.NInteractionTerms(Ord1,Nc)

  N,c  = size(MV1)
  Nby2 = Int(N/2)

  h    = zeros(T,N)

  topi = 1:Nby2
  boti = Nby2+1:N

  if (Ord1 + 1 == Ord0) && Ord0>=2

    # The Laplacian nonlinearity {δ4(∇A); δ4'∇(A*)}
    for i in 1:Nt1
      ind1 = CenterManifold.GetPolynomialIndices(i,Ord1,Nc) .+ 1

      ind_total = [ind1[:]; mode]
      sort!(ind_total)
      if (ind_total == ind0)    # Matching polynomials
        # println("Index Match: $ind1, $mode: $ind0")
        h[topi] .= h[topi] .+ Lap*MV1[topi,i]
      end

      ind_total = [ind1[:]; modeC]
      sort!(ind_total)
      if (ind_total == ind0)    # Matching polynomials
        # println("Index Match: $ind1, $modeC: $ind0")
        h[boti] .= h[boti] .+ LapC*MV1[boti,i]
      end            
    end           # i
  end             # if (Ord1 + 1 == Ord0) && Ord0>=2

  return h
end
#---------------------------------------------------------------------- 
function BuildLowOrder(n0::Int,Ord::Int,Nc::Int,Y::AbstractMatrix{T},G::AbstractMatrix{T},OrdY::Int,OrdG::Int) where {T<:Number}

  r,c       = size(Y)

  h_low     = zeros(T,r)

  if (OrdY-1+OrdG != Ord)
    return h_low
  end

  NtY       = CenterManifold.NInteractionTerms(OrdY,Nc)
  NtG       = CenterManifold.NInteractionTerms(OrdG,Nc)

  for α in 1:Nc         # Derivative Variable
    for ii in 1:OrdY    # Polynomial order
      for i in 1:NtY    # Go through all terms that undergo derivation

        ind = CenterManifold.GetPolynomialIndices(i,OrdY,Nc)

        # If ∂/∂x_α is non-zero. 
        if (ind[ii] == (α-1))
          ind2 = [ind[1:ii-1]; ind[ii+1:end]]
          
          # Multiplication by G_α,j
          for j in 1:NtG
            ind3        = CenterManifold.GetPolynomialIndices(i,OrdG,Nc)
            ind4        = [ind2[:]; ind3] # Indicies of Y_i * G_α,j

            k           = CenterManifold.GetIndexNumber(ind4,Nc)
            if k == n0
              h_low    .= h_low .+ G[α,j]*Y[:,i]
            end  
          end     # j
        end       # if (ind[ii] == α-1)
      end         # i
    end           # ii
  end             # α

  return h_low
end  
#----------------------------------------------------------------------
function GL2AsymptoticCM_Ord2(L::AbstractMatrix{T2},LC::AbstractMatrix{T2},B::AbstractVector{T3},δ::AbstractVector{T2},Khat::AbstractMatrix{T1},Vext::AbstractMatrix{T2},Wext::AbstractMatrix{T2},n::Int,p::Int,h::Int,lbc::Bool,rbc::Bool,Lap::AbstractMatrix{T2},LapC::AbstractMatrix{T2},ifnormal::Bool) where {T1,T2,T3<:Number}

  # n         - Original Critical space size (Needed to determine Laplace Modes)
  # p         - Parameter Modes
  # h         - Forcing Modes
  Ord         = 2
  Nby2        = size(L,2)
  N           = Nby2*2
  m           = size(Khat,2)
  ind1        = 1:Nby2
  ind2        = Nby2+1:N
  restricted  = ~ifnormal

  # No of terms at this order.
  Nt          = CenterManifold.NInteractionTerms(Ord,m)
  # Asymptotic CM Vectors.
  Y           = zeros(T2,N,Nt)
  # Asymptotic Forcing Vectors.
  F           = zeros(T2,N,Nt)
  # Lower Triangular Linear System Matrix at this order for the asymptotic vectors Y
  SysMat      = CenterManifold.BuildAsympSystem(Ord,Khat)
  # Reduced Matrix terms.
  G           = zeros(T1,m,Nt)


  # Stepper-Arnoldi (for Extended Operator Calculation)
  #-------------------------------------------------- 
  ifadjoint         = false
  ifoptimal         = false
  ifverbose         = false
  verbosestep       = 500
  nsteps            = 500
  dt                = 1.0e-4
  SInp              = StepperArnoldi.StepperInput(ifadjoint,ifoptimal,ifverbose,verbosestep,nsteps,dt)
  
  ifarnoldi         = true 
  ifverbose         = false
  ifeigshift        = true
  vlen              = Nby2
  nev               = 1
  ekryl             = 15  
  lkryl             = nev + ekryl 
  eigshift          = 0.0 + 0.0im
  ngs               = 2
  bsize             = 1
  outer_iterations  = 100
  tol               = 1.0e-12
  AInp              = StepperArnoldi.ArnoldiInput(ifarnoldi,ifverbose,ifeigshift,vlen,nev,ekryl,lkryl,eigshift,ngs,bsize,outer_iterations,tol)
  #-------------------------------------------------- 

  # Build the forcing
  for i in 1:Nt
    ω       = SysMat[i,i]

    println("Solving for $(i), ω=$(ω)")

    h_asymp = zeros(T2,N)
    if p>0
      lap_mode  = n+1
      lap_cmode = n+2
      Ord1      = 1
      h_lapl    = GLLapNLTerm(i,Ord,m,Vext,Ord1,lap_mode,lap_cmode,BiLapg,BiLapCg)
      h_asymp  .= h_asymp .+ h_lapl 
    end

    Ord1        = 1
    Ord2        = 1
    Ord3        = 1
    δ5          = δ[5]
    h_NL        = GLStdAsympNLTerm(i,Ord,m,Vext,Vext,Vext,Ord1,Ord2,Ord3,δ5)
    h_asymp    .= h_asymp .+ h_NL

    h_offd      = zeros(T2,N)
    for j in 1:(i-1)
      h_offd .= h_offd .+ SysMat[i,j]*Y[:,j]
    end
    h_asymp  .= h_asymp .- h_offd
    @views SEM1D.SEM_SetBC!(h_asymp[ind1],lbc,rbc)
    @views SEM1D.SEM_SetBC!(h_asymp[ind2],lbc,rbc)
    for j in eachindex(h_asymp)
      F[j,i]  = h_asymp[j]
    end

    V1        = view(Vext,ind1,:)
    V2        = view(Vext,ind2,:)
    W1        = view(Wext,ind1,:)
    W2        = view(Wext,ind2,:)
    h_asymp1  = view(h_asymp,ind1)
    h_asymp2  = view(h_asymp,ind2)

    λ_khat    = diag(Khat)
    AInp.eigshift = ω
    if (restricted)
      CM2_EM1   = StepperArnoldi.ExtendTangentSpaceRestricted(L,B,λ_khat,V1,W1,h_asymp1,ω,AInp,SInp,lbc,rbc)
      CM2_EM2   = StepperArnoldi.ExtendTangentSpaceRestricted(LC,B,λ_khat,V2,W2,h_asymp2,ω,AInp,SInp,lbc,rbc)
    else
      CM2_EM1   = StepperArnoldi.ExtendTangentSpace(L,B,λ_khat,V1,W1,h_asymp1,ω,AInp,SInp,lbc,rbc)
      CM2_EM2   = StepperArnoldi.ExtendTangentSpace(LC,B,λ_khat,V2,W2,h_asymp2,ω,AInp,SInp,lbc,rbc)
    end

    Y[ind1,i] = copy(CM2_EM1.ve)
    Y[ind2,i] = copy(CM2_EM2.ve)

    for j in 1:m
      G[j,i]  = CM2_EM1.Γ[j] + CM2_EM2.Γ[j]
    end     # j

  end       # i

  return SysMat,G,Y,F
end  
#---------------------------------------------------------------------- 
function GL2AsymptoticCM_Ord3(L::AbstractMatrix{T2},LC::AbstractMatrix{T2},B::AbstractVector{T3},δ::AbstractVector{T2},Khat::AbstractMatrix{T1},Vext::AbstractMatrix{T2},Wext::AbstractMatrix{T2},G2::AbstractMatrix{T1},Y2::AbstractMatrix{T2},n::Int,p::Int,h::Int,lbc::Bool,rbc::Bool,LAP::AbstractMatrix{T2},LAPC::AbstractMatrix{T2},irnormal::Bool) where {T1,T2,T3<:Number}

  # n         - Original Critical space size (Needed to determine Laplace Modes)
  # p         - Parameter Modes
  # h         - Forcing Modes
  Ord         = 3
  Nby2        = size(L,2)
  N           = Nby2*2
  m           = size(Khat,2)
  ind1        = 1:Nby2
  ind2        = Nby2+1:N
  restricted  = ~ifnormal

  # No of terms at this order.
  Nt          = CenterManifold.NInteractionTerms(Ord,m)
  # Asymptotic CM Vectors.
  Y           = zeros(T2,N,Nt)
  # Asymptotic Forcing Vectors.
  F           = zeros(T2,N,Nt)
  # Lower Triangular Linear System Matrix at this order for the asymptotic vectors Y
  SysMat      = CenterManifold.BuildAsympSystem(Ord,Khat)
  # Reduced Matrix terms.
  G           = zeros(T1,m,Nt)



  # Stepper-Arnoldi (for Extended Operator Calculation)
  #-------------------------------------------------- 
  ifadjoint         = false
  ifoptimal         = false
  ifverbose         = false
  verbosestep       = 500
  nsteps            = 500
  dt                = 1.0e-4
  SInp              = StepperArnoldi.StepperInput(ifadjoint,ifoptimal,ifverbose,verbosestep,nsteps,dt)
  
  ifarnoldi         = true 
  ifverbose         = false
  ifeigshift        = true
  vlen              = Nby2
  nev               = 1
  ekryl             = 15  
  lkryl             = nev + ekryl 
  eigshift          = 0.0 + 0.0im
  ngs               = 2
  bsize             = 1
  outer_iterations  = 100
  tol               = 1.0e-12
  AInp              = StepperArnoldi.ArnoldiInput(ifarnoldi,ifverbose,ifeigshift,vlen,nev,ekryl,lkryl,eigshift,ngs,bsize,outer_iterations,tol)
  #-------------------------------------------------- 

  # Build the forcing
  for i in 1:Nt
    ω       = SysMat[i,i]

    println("Solving for $(i), ω=$(ω)")

    h_asymp = zeros(T2,N)
    if p>0
      lap_mode  = n+1
      lap_cmode = n+2
      OrdY      = 2
      h_lapl    = GLLapNLTerm(i,Ord,m,Y2,OrdY,lap_mode,lap_cmode,LAP,LAPC)
      h_asymp  .= h_asymp .+ h_lapl 
    end

    Ord1        = 1
    Ord2        = 1
    Ord3        = 1
    δ5          = δ[5]
    h_NL        = GLStdAsympNLTerm(i,Ord,m,Vext,Vext,Vext,Ord1,Ord2,Ord3,δ5)
    h_asymp    .= h_asymp .+ h_NL

    h_offd      = zeros(T2,N)
    for j in 1:(i-1)
      h_offd .= h_offd .+ SysMat[i,j]*Y[:,j]
    end
    h_asymp  .= h_asymp .- h_offd

    OrdY      = 2
    OrdG      = 2
    h_low     = BuildLowOrder(i,Ord,m,Y2,G2,OrdY,OrdG)
    h_asymp  .= h_asymp .- h_low

    @views SEM1D.SEM_SetBC!(h_asymp[ind1],lbc,rbc)
    @views SEM1D.SEM_SetBC!(h_asymp[ind2],lbc,rbc)
    for j in eachindex(h_asymp)
      F[j,i]  = h_asymp[j]
    end

    V1        = view(Vext,ind1,:)
    V2        = view(Vext,ind2,:)
    W1        = view(Wext,ind1,:)
    W2        = view(Wext,ind2,:)
    h_asymp1  = view(h_asymp,ind1)
    h_asymp2  = view(h_asymp,ind2)

    λ_khat    = diag(Khat)
    AInp.eigshift = ω
    if (restricted)
      CM2_EM1   = StepperArnoldi.ExtendTangentSpaceRestricted(L,B,λ_khat,V1,W1,h_asymp1,ω,AInp,SInp,lbc,rbc)
      CM2_EM2   = StepperArnoldi.ExtendTangentSpaceRestricted(LC,B,λ_khat,V2,W2,h_asymp2,ω,AInp,SInp,lbc,rbc)
    else
      CM2_EM1   = StepperArnoldi.ExtendTangentSpace(L,B,λ_khat,V1,W1,h_asymp1,ω,AInp,SInp,lbc,rbc)
      CM2_EM2   = StepperArnoldi.ExtendTangentSpace(LC,B,λ_khat,V2,W2,h_asymp2,ω,AInp,SInp,lbc,rbc)
    end

    Y[ind1,i] = copy(CM2_EM1.ve)
    Y[ind2,i] = copy(CM2_EM2.ve)

    for j in 1:m
      G[j,i]  = CM2_EM1.Γ[j] + CM2_EM2.Γ[j]
    end     # j

  end       # i

  return SysMat,G,Y,F
end  
#---------------------------------------------------------------------- 
resmodeplot = true

SM2,G2,Y2,F2 = GL2AsymptoticCM_Ord2(OPg,OPCg,Bg,δ,Khat,Vext,Wext,n,p,h,Inp.lbc,Inp.rbc,BiLapg,BiLapCg,ifnormal);
if (resmodeplot)
  h4    = figure(num=4,figsize=Grh.figsz2)
  ax4   = gca()
  ax4.cla()
  plno  = -1 
end

for i in 1:CenterManifold.NInteractionTerms(2,m)
  vnorm = norm(Y2[ind1,i])
  # Plot Mode
  if (resmodeplot) && vnorm > 0.0
    j = i
    global plno = plno + 1

    ind = CenterManifold.GetPolynomialIndices(i,2,m) .+ 1
    ijk = "_{$(ind[1])"
    for k in 2:2
      ijk = ijk*",$(ind[k])"
    end
    ijk = ijk*"}"
    # leg = L"\mathfrak{R}(y%$ijk)"
    ax4.plot(xg,real.(Y2[ind1,i]),linewidth=2,linestyle="-", color=cm(plno),label=L"\mathfrak{R}(y_{%$ijk})")
    ax4.plot(xg,imag.(Y2[ind1,i]),linewidth=2,linestyle="--",color=cm(plno),label=L"\mathfrak{Im}(y_{%$ijk})")
  end
end  
if (resmodeplot)
  ax4.legend(ncols=4,fontsize=Grh.lgfs)
end

SM3,G3,Y3,F3 = GL2AsymptoticCM_Ord3(OPg,OPCg,Bg,δ,Khat,Vext,Wext,G2,Y2,n,p,h,Inp.lbc,Inp.rbc,BiLapg,BiLapCg,ifnormal);

if (resmodeplot)
  h5    = figure(num=5,figsize=Grh.figsz2)
  ax5   = gca()
  ax5.cla()
  plno  = -1 
end

for i in 1:CenterManifold.NInteractionTerms(3,m)
  vnorm = norm(Y3[ind1,i])
  # Plot Mode
  if (resmodeplot) && vnorm > 0.0
 
    ind = CenterManifold.GetPolynomialIndices(i,3,m) .+ 1
    ijk = "_{$(ind[1])"
    for k in 2:3
      ijk = ijk*",$(ind[k])"
    end
    ijk = ijk*"}"
   
    j = i
    global plno = plno + 1 
    ax5.plot(xg,real.(Y3[ind1,i]),linewidth=2,linestyle="-", color=cm(plno),label=L"\mathfrak{R}(y_{%$ijk})")
    ax5.plot(xg,imag.(Y3[ind1,i]),linewidth=2,linestyle="--",color=cm(plno),label=L"\mathfrak{Im}(y_{%$ijk})")
  end
end  
if (resmodeplot)
  ax5.legend(ncols=4,fontsize=Grh.lgfs,loc="upper right")
end

println("Asymptotic System (Normal Form) Done.")
















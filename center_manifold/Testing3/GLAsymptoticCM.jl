# Calculate the Center-Manifold Asymptotic vectors
include("../Module_CenterManifold/CenterManifold.jl")
#---------------------------------------------------------------------- 
include("GLAsymptoticNLTerms.jl")
#----------------------------------------------------------------------
function GL2AsymptoticCM_Ord2(L::AbstractMatrix{T2},LC::AbstractMatrix{T2},B::AbstractVector{T3},δ::AbstractVector{T2},Khat::AbstractMatrix{T1},Vext::AbstractMatrix{T2},Wext::AbstractMatrix{T2},σ::AbstractVector{T1},PertModes::AbstractVector{Int},nsys::Int,p::Int,h::Int,lbc::Bool,rbc::Bool,Lap::AbstractMatrix{T2},LapC::AbstractMatrix{T2},ifnormal::Bool) where {T1,T2,T3<:Number}

  # nsys      - Original (Perturbed) Critical space size (Needed to determine Laplace Modes and Perturbation Modes)
  # p         - Parameter Modes
  # h         - Forcing Modes
  Ord         = 2
  Nby2        = size(L,2)
  N           = Nby2*2
  m           = size(Khat,2)
  ind1        = 1:Nby2
  ind2        = Nby2+1:N
  restricted  = ~ifnormal

  PM          = PertModes[PertModes .> 0]
  npmodes     = length(PM)
  # println(npmodes)
  # println(PertModes)

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
  # Mass Matrix (Vector)
  B2          = [B; B]

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

  AsEM1             = Vector{StepperArnoldi.ExtendedMode}(undef,Nt)
  AsEM2             = Vector{StepperArnoldi.ExtendedMode}(undef,Nt)
  Z                 = zeros(T2,Nt,m)
  Γ                 = zeros(T2,m,Nt)
  We                = zeros(T2,N,Nt)

  # Build the forcing
  for i in 1:Nt
    ω       = SysMat[i,i]

    println("Solving for $(i), ω=$(ω)")

    h_asymp = zeros(T2,N)
    if npmodes>0
      Ord1      = 1
      h_mpert   = GLModePertTerm(i,Ord,Vext,Ord1,Vext,Wext,B2,PertModes)
      # if norm(h_mpert)>1.0e-12
      #   println("Pertmode for i=$i, |h_mpert| = $(norm(h_mpert))")
      # end
      h_asymp  .= h_asymp .+ h_mpert
    end

    if p>0
      lap_mode  = nsys+npmodes + 1
      lap_cmode = nsys+npmodes + 2
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

    AInp.ifeigshift = true
    AInp.eigshift   = ω

    if (npmodes>0)
      AsEM1[i]    = StepperArnoldi.REPTangentSpace(L, B,λ_khat,σ,V1,W1,h_asymp1,ω,restricted,AInp,SInp,lbc,rbc)
      AsEM2[i]    = StepperArnoldi.REPTangentSpace(LC,B,λ_khat,σ,V2,W2,h_asymp2,ω,restricted,AInp,SInp,lbc,rbc)
    else
      AsEM1[i]    = StepperArnoldi.ExtendTangentSpace(L, B,λ_khat,V1,W1,h_asymp1,ω,restricted,AInp,SInp,lbc,rbc)
      AsEM2[i]    = StepperArnoldi.ExtendTangentSpace(LC,B,λ_khat,V2,W2,h_asymp2,ω,restricted,AInp,SInp,lbc,rbc)
    end

    Y[ind1,i]     = copy(AsEM1[i].ve)
    Y[ind2,i]     = copy(AsEM2[i].ve)
    Z[i,:]        = AsEM1[i].z .+ AsEM2[i].z
    Γ[:,i]        = AsEM1[i].Γ .+ AsEM2[i].Γ

    We[ind1,i]    = copy(AsEM1[i].we)
    We[ind2,i]    = copy(AsEM2[i].we)

    for j in 1:m
      G[j,i]  = AsEM1[i].Γ[j] + AsEM2[i].Γ[j]
    end     # j

  end       # i

  λ_khat  = diag(Khat)
  CMmodes = StepperArnoldi.ExtendedModes(λ_khat,F,Γ,Y,We,Z)

  return SysMat,G,Y,F,CMmodes
end  
#---------------------------------------------------------------------- 
function GL2AsymptoticCM_Ord3(L::AbstractMatrix{T2},LC::AbstractMatrix{T2},B::AbstractVector{T3},δ::AbstractVector{T2},Khat::AbstractMatrix{T1},Vext::AbstractMatrix{T2},Wext::AbstractMatrix{T2},σ::AbstractVector{T1},PertModes::AbstractVector{Int},nsys::Int,p::Int,h::Int,G2::AbstractMatrix{T1},Y2::AbstractMatrix{T2},lbc::Bool,rbc::Bool,LAP::AbstractMatrix{T2},LAPC::AbstractMatrix{T2},ifnormal::Bool) where {T1,T2,T3<:Number}

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

  PM          = PertModes[PertModes .> 0]
  npmodes     = length(PM)

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
  # Mass Matrix (Vector)
  B2          = [B; B]


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

  AsEM1             = Vector{StepperArnoldi.ExtendedMode}(undef,Nt)
  AsEM2             = Vector{StepperArnoldi.ExtendedMode}(undef,Nt)
  Z                 = zeros(T2,Nt,m)
  Γ                 = zeros(T2,m,Nt)
  We                = zeros(T2,N,Nt)

  # Build the forcing
  for i in 1:Nt
    ω       = SysMat[i,i]

    println("Solving for $(i), ω=$(ω)")

    h_asymp = zeros(T2,N)
    if npmodes>0
      Ord1      = 2
      h_mpert   = GLModePertTerm(i,Ord,Y2,Ord1,Vext,Wext,B2,PertModes)
      # if norm(h_mpert)>1.0e-12
      #   println("Pertmode for i=$i, |h_mpert| = $(norm(h_mpert))")
      # end
      h_asymp  .= h_asymp .+ h_mpert
    end

    if p>0
      lap_mode  = nsys+npmodes + 1
      lap_cmode = nsys+npmodes + 2
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

    λ_khat          = diag(Khat)
    AInp.ifeigshift = true
    AInp.eigshift   = ω

    if (npmodes>0)
      AsEM1[i]    = StepperArnoldi.REPTangentSpace(L, B,λ_khat,σ,V1,W1,h_asymp1,ω,restricted,AInp,SInp,lbc,rbc)
      AsEM2[i]    = StepperArnoldi.REPTangentSpace(LC,B,λ_khat,σ,V2,W2,h_asymp2,ω,restricted,AInp,SInp,lbc,rbc)
    else
      AsEM1[i]    = StepperArnoldi.ExtendTangentSpace(L, B,λ_khat,V1,W1,h_asymp1,ω,restricted,AInp,SInp,lbc,rbc)
      AsEM2[i]    = StepperArnoldi.ExtendTangentSpace(LC,B,λ_khat,V2,W2,h_asymp2,ω,restricted,AInp,SInp,lbc,rbc)
    end

    Y[ind1,i]     = copy(AsEM1[i].ve)
    Y[ind2,i]     = copy(AsEM2[i].ve)
    Z[i,:]        = AsEM1[i].z .+ AsEM2[i].z
    Γ[:,i]        = AsEM1[i].Γ .+ AsEM2[i].Γ

    We[ind1,i]    = copy(AsEM1[i].we)
    We[ind2,i]    = copy(AsEM2[i].we)

    for j in 1:m
      G[j,i]  = AsEM1[i].Γ[j] + AsEM2[i].Γ[j]
    end     # j

  end       # i

  λ_khat  = diag(Khat)
  CMmodes = StepperArnoldi.ExtendedModes(λ_khat,F,Γ,Y,We,Z)

  return SysMat,G,Y,F,CMmodes
end  
#---------------------------------------------------------------------- 
function GL2AsymptoticCM_Ord4(L::AbstractMatrix{T2},LC::AbstractMatrix{T2},B::AbstractVector{T3},δ::AbstractVector{T2},Khat::AbstractMatrix{T1},Vext::AbstractMatrix{T2},Wext::AbstractMatrix{T2},σ::AbstractVector{T1},PertModes::AbstractVector{Int},nsys::Int,p::Int,h::Int,G2::AbstractMatrix{T1},Y2::AbstractMatrix{T2},G3::AbstractMatrix{T1},Y3::AbstractMatrix{T2},lbc::Bool,rbc::Bool,LAP::AbstractMatrix{T2},LAPC::AbstractMatrix{T2},ifnormal::Bool) where {T1,T2,T3<:Number}

  # nsys      - Original Critical space size (Needed to determine Laplace Modes)
  # p         - Parameter Modes
  # h         - Forcing Modes
  Ord         = 4
  Nby2        = size(L,2)
  N           = Nby2*2
  m           = size(Khat,2)
  ind1        = 1:Nby2
  ind2        = Nby2+1:N
  restricted  = ~ifnormal

  PM          = PertModes[PertModes .> 0]
  npmodes     = length(PM)

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
  # Mass Matrix (Vector)
  B2          = [B; B]


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

  AsEM1             = Vector{StepperArnoldi.ExtendedMode}(undef,Nt)
  AsEM2             = Vector{StepperArnoldi.ExtendedMode}(undef,Nt)
  Z                 = zeros(T2,Nt,m)
  Γ                 = zeros(T2,m,Nt)
  We                = zeros(T2,N,Nt)

  # Build the forcing
  for i in 1:Nt
    ω       = SysMat[i,i]

    println("Solving for $(i), ω=$(ω)")
    h_asymp = zeros(T2,N)

    # Mode Perturbation Terms
    #-------------------------------------------------- 
    if npmodes>0
      Ord1      = 3
      h_mpert   = GLModePertTerm(i,Ord,Y3,Ord1,Vext,Wext,B2,PertModes)
      h_asymp  .= h_asymp .+ h_mpert
    end

    # Parameter Perturbation Terms
    #--------------------------------------------------  
    if p>0
      lap_mode  = nsys+npmodes + 1
      lap_cmode = nsys+npmodes + 2
      OrdY      = 3
      h_lapl    = GLLapNLTerm(i,Ord,m,Y3,OrdY,lap_mode,lap_cmode,LAP,LAPC)
      h_asymp  .= h_asymp .+ h_lapl 
    end

    # Standard Non-Linearity
    #--------------------------------------------------  
    Ord1        = 1
    Ord2        = 1
    Ord3        = 2
    δ5          = δ[5]
    h_NL        = GLStdAsympNLTerm(i,Ord,m,Vext,Vext,Y2,Ord1,Ord2,Ord3,δ5)
    h_asymp    .= h_asymp .+ h_NL

    # Off-Diagonal Terms
    #--------------------------------------------------  
    h_offd      = zeros(T2,N)
    for j in 1:(i-1)
      h_offd .= h_offd .+ SysMat[i,j]*Y[:,j]
    end
    h_asymp  .= h_asymp .- h_offd

    # Low Order Terms
    #-------------------------------------------------- 
    OrdY      = 2
    OrdG      = 3
    h_low     = BuildLowOrder(i,Ord,m,Y2,G3,OrdY,OrdG)
    h_asymp  .= h_asymp .- h_low

    OrdY      = 3
    OrdG      = 2
    h_low     = BuildLowOrder(i,Ord,m,Y3,G2,OrdY,OrdG)
    h_asymp  .= h_asymp .- h_low
    #--------------------------------------------------   


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

    λ_khat          = diag(Khat)
    AInp.ifeigshift = true
    AInp.eigshift   = ω

    if (npmodes>0)
      AsEM1[i]    = StepperArnoldi.REPTangentSpace(L, B,λ_khat,σ,V1,W1,h_asymp1,ω,restricted,AInp,SInp,lbc,rbc)
      AsEM2[i]    = StepperArnoldi.REPTangentSpace(LC,B,λ_khat,σ,V2,W2,h_asymp2,ω,restricted,AInp,SInp,lbc,rbc)
    else
      AsEM1[i]    = StepperArnoldi.ExtendTangentSpace(L, B,λ_khat,V1,W1,h_asymp1,ω,restricted,AInp,SInp,lbc,rbc)
      AsEM2[i]    = StepperArnoldi.ExtendTangentSpace(LC,B,λ_khat,V2,W2,h_asymp2,ω,restricted,AInp,SInp,lbc,rbc)
    end

    Y[ind1,i]     = copy(AsEM1[i].ve)
    Y[ind2,i]     = copy(AsEM2[i].ve)
    Z[i,:]        = AsEM1[i].z .+ AsEM2[i].z
    Γ[:,i]        = AsEM1[i].Γ .+ AsEM2[i].Γ

    We[ind1,i]    = copy(AsEM1[i].we)
    We[ind2,i]    = copy(AsEM2[i].we)

    for j in 1:m
      G[j,i]  = AsEM1[i].Γ[j] + AsEM2[i].Γ[j]
    end     # j

  end       # i

  λ_khat  = diag(Khat)
  CMmodes = StepperArnoldi.ExtendedModes(λ_khat,F,Γ,Y,We,Z)

  return SysMat,G,Y,F,CMmodes
end  
#---------------------------------------------------------------------- 
function GL2AsymptoticCM_Ord5(L::AbstractMatrix{T2},LC::AbstractMatrix{T2},B::AbstractVector{T3},δ::AbstractVector{T2},Khat::AbstractMatrix{T1},Vext::AbstractMatrix{T2},Wext::AbstractMatrix{T2},σ::AbstractVector{T1},PertModes::AbstractVector{Int},nsys::Int,p::Int,h::Int,G2::AbstractMatrix{T1},Y2::AbstractMatrix{T2},G3::AbstractMatrix{T1},Y3::AbstractMatrix{T2},G4::AbstractMatrix{T1},Y4::AbstractMatrix{T2},lbc::Bool,rbc::Bool,LAP::AbstractMatrix{T2},LAPC::AbstractMatrix{T2},ifnormal::Bool) where {T1,T2,T3<:Number}

  # nsys      - Original Critical space size (Needed to determine Laplace Modes)
  # p         - Parameter Modes
  # h         - Forcing Modes
  Ord         = 5
  Nby2        = size(L,2)
  N           = Nby2*2
  m           = size(Khat,2)
  ind1        = 1:Nby2
  ind2        = Nby2+1:N
  restricted  = ~ifnormal

  PM          = PertModes[PertModes .> 0]
  npmodes     = length(PM)

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
  # Mass Matrix (Vector)
  B2          = [B; B]


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

  AsEM1             = Vector{StepperArnoldi.ExtendedMode}(undef,Nt)
  AsEM2             = Vector{StepperArnoldi.ExtendedMode}(undef,Nt)
  Z                 = zeros(T2,Nt,m)
  Γ                 = zeros(T2,m,Nt)
  We                = zeros(T2,N,Nt)

  # Build the forcing
  for i in 1:Nt
    ω       = SysMat[i,i]

    println("Solving for $(i), ω=$(ω)")

    # Mode Perturbation Terms
    #-------------------------------------------------- 
    h_asymp = zeros(T2,N)
    if npmodes>0
      Ord1      = 4
      h_mpert   = GLModePertTerm(i,Ord,Y4,Ord1,Vext,Wext,B2,PertModes)
      h_asymp  .= h_asymp .+ h_mpert
    end

    # Parameter Perturbation Terms
    #--------------------------------------------------  
    if p>0
      lap_mode  = nsys+npmodes + 1
      lap_cmode = nsys+npmodes + 2
      OrdY      = 4
      h_lapl    = GLLapNLTerm(i,Ord,m,Y4,OrdY,lap_mode,lap_cmode,LAP,LAPC)
      h_asymp  .= h_asymp .+ h_lapl 
    end

    # Standard Non-Linearity
    #--------------------------------------------------  
    Ord1        = 1
    Ord2        = 1
    Ord3        = 3
    δ5          = δ[5]
    h_NL        = GLStdAsympNLTerm(i,Ord,m,Vext,Vext,Y3,Ord1,Ord2,Ord3,δ5)
    h_asymp    .= h_asymp .+ h_NL

    Ord1        = 1
    Ord2        = 2
    Ord3        = 2
    δ5          = δ[5]
    h_NL        = GLStdAsympNLTerm(i,Ord,m,Vext,Y2,Y2,Ord1,Ord2,Ord3,δ5)
    h_asymp    .= h_asymp .+ h_NL


    # Off-Diagonal Terms
    #--------------------------------------------------  
    h_offd      = zeros(T2,N)
    for j in 1:(i-1)
      h_offd .= h_offd .+ SysMat[i,j]*Y[:,j]
    end
    h_asymp  .= h_asymp .- h_offd

    # Low Order Terms
    #-------------------------------------------------- 
    OrdY      = 2
    OrdG      = 4
    h_low     = BuildLowOrder(i,Ord,m,Y2,G4,OrdY,OrdG)
    h_asymp  .= h_asymp .- h_low

    OrdY      = 3
    OrdG      = 3
    h_low     = BuildLowOrder(i,Ord,m,Y3,G3,OrdY,OrdG)
    h_asymp  .= h_asymp .- h_low

    OrdY      = 4
    OrdG      = 2
    h_low     = BuildLowOrder(i,Ord,m,Y4,G2,OrdY,OrdG)
    h_asymp  .= h_asymp .- h_low
    #--------------------------------------------------   

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

    λ_khat          = diag(Khat)
    AInp.ifeigshift = true
    AInp.eigshift   = ω

    if (npmodes>0)
      AsEM1[i]    = StepperArnoldi.REPTangentSpace(L, B,λ_khat,σ,V1,W1,h_asymp1,ω,restricted,AInp,SInp,lbc,rbc)
      AsEM2[i]    = StepperArnoldi.REPTangentSpace(LC,B,λ_khat,σ,V2,W2,h_asymp2,ω,restricted,AInp,SInp,lbc,rbc)
    else
      AsEM1[i]    = StepperArnoldi.ExtendTangentSpace(L, B,λ_khat,V1,W1,h_asymp1,ω,restricted,AInp,SInp,lbc,rbc)
      AsEM2[i]    = StepperArnoldi.ExtendTangentSpace(LC,B,λ_khat,V2,W2,h_asymp2,ω,restricted,AInp,SInp,lbc,rbc)
    end

    Y[ind1,i]     = copy(AsEM1[i].ve)
    Y[ind2,i]     = copy(AsEM2[i].ve)
    Z[i,:]        = AsEM1[i].z .+ AsEM2[i].z
    Γ[:,i]        = AsEM1[i].Γ .+ AsEM2[i].Γ

    We[ind1,i]    = copy(AsEM1[i].we)
    We[ind2,i]    = copy(AsEM2[i].we)

    for j in 1:m
      G[j,i]  = AsEM1[i].Γ[j] + AsEM2[i].Γ[j]
    end     # j

  end       # i

  λ_khat  = diag(Khat)
  CMmodes = StepperArnoldi.ExtendedModes(λ_khat,F,Γ,Y,We,Z)

  return SysMat,G,Y,F,CMmodes
end  
#---------------------------------------------------------------------- 



#---------------------------------------------------------------------- 
println("Asymptotic System (Normal Form) Done.")
















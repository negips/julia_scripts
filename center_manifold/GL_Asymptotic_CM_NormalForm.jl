# Testing the module

# include("Module_SEM1D/SEM1D.jl")
# #using .SEM1D
# 
# include("Module_StepperArnoldi/StepperArnoldi.jl")
# #using .StepperArnoldi
# 
# include("Module_CenterManifold/CenterManifold.jl")
# #using .CenterManifold

#---------------------------------------------------------------------- 
function GLStdAsympNLTerm(n0::Int,Ord0::Int,nc::Int,MV1::AbstractMatrix{T},MV2::AbstractMatrix{T},MV3::AbstractMatrix{T},Ord1::Int,Ord2::Int,Ord3::Int) where {T}

  ind0 = CenterManifold.GetPolynomialIndices(n0,Ord0,nc) .+ 1

  Nt1  = CenterManifold.NInteractionTerms(Ord1,nc)
  Nt2  = CenterManifold.NInteractionTerms(Ord2,nc)
  Nt3  = CenterManifold.NInteractionTerms(Ord3,nc)

  N,c  = size(MV1)
  Nby2 = Int(N/2)

  h    = zeros(T,N)

  topi = 1:Nby2
  boti = Nby2+1:N

  if (Ord1 + Ord2 + Ord3 == Ord0) && Ord0>=3

    # The standard nonlinearity δ5*|A*A|A and δ5'|A*A|A* only starts at third-order.
    for i in 1:Nt1
      ind1 = CenterManifold.GetPolynomialIndices(i,Ord1,nc) .+ 1
      for j in 1:Nt2
        ind2 = CenterManifold.GetPolynomialIndices(j,Ord2,nc) .+ 1
        for k in 1:Nt3
          ind3 = CenterManifold.GetPolynomialIndices(k,Ord3,nc) .+ 1

          ind_total = [ind1[:];ind2[:];ind3[:]]

          sort!(ind_total)

          if (ind_total == ind0)    # Matching polynomials
            # println("Index Match: $ind1, $ind2, $ind3: $ind0")
            h[topi] .= h[topi] .+ MV1[boti,i].*MV2[topi,j].*MV3[topi,k]
            h[boti] .= h[boti] .+ MV1[topi,i].*MV2[boti,j].*MV3[boti,k]
          end            
        end       # k
      end         # j
    end           # i
  end             # if (Ord1 + Ord2 + Ord3 == Ord0) && Ord0>=3

  return h
end
#---------------------------------------------------------------------- 
function GLLapNLTerm(n0::Int,Ord0::Int,nc::Int,MV1::AbstractMatrix{T},Ord1::Int,mode::Int,modeC::Int,Lap::AbstractMatrix{S},LapC::AbstractMatrix{S}) where {T,S}

  # mode    - Mode No of δ4
  # modeC   - Mode No of δ4' (conjugate)

  ind0 = CenterManifold.GetPolynomialIndices(n0,Ord0,nc) .+ 1

  Nt1  = CenterManifold.NInteractionTerms(Ord1,nc)

  N,c  = size(MV1)
  Nby2 = Int(N/2)

  h    = zeros(T,N)

  topi = 1:Nby2
  boti = Nby2+1:N

  if (Ord1 + 1 == Ord0) && Ord0>=2

    # The Laplacian nonlinearity {δ4(∇A); δ4'∇(A*)}
    for i in 1:Nt1
      ind1 = CenterManifold.GetPolynomialIndices(i,Ord1,nc) .+ 1

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
function BuildAsympSystem(Ord::Int,Nc::Int,Khat::AbstractMatrix{T}) where {T<:Number}

  Nt        = CenterManifold.NInteractionTerms(Ord,Nc)
  SysM      = zeros(T,Nt,Nt)
  for α in 1:Nc         # Derivative Variable
    for ii in 1:Ord     # Polynomial order
      for i in 1:Nt     # Go through all terms that undergo derivation

        ind = CenterManifold.GetPolynomialIndices(i,Ord,Nc)

        # If ∂/∂x_α is non-zero. 
        if (ind[ii] == (α-1))
          ind2 = [ind[1:ii-1]; ind[ii+1:end]]

          # Multiplication by K_α,j
          for j in 1:Nc
            ind3        = [ind2[:]; (j-1)]
            k           = CenterManifold.GetIndexNumber(ind3,Nc)
            SysM[k,i]   = SysM[k,i] + Khat[α,j]
          end     # j
        end       # if (ind[ii] == α-1)
      end         # i
    end           # ii
  end             # α

  return SysM
end  
#----------------------------------------------------------------------

# Second Order Asymptotic terms
Ord         = 2
Nt          = CenterManifold.NInteractionTerms(Ord,m)
Y_O2        = zeros(ComplexF64,N,Nt)
SysMat_O2   = BuildAsympSystem(Ord,m,Khat)

# Reduced Matrix terms
G_O2        = zeros(ComplexF64,m,Nt)

for i in 1:Nt
  ind       = CenterManifold.GetPolynomialIndices(i,Ord,m)

  println("Solving for $(ind .+ 1)")

  h_asymp   = zeros(ComplexF64,N)
  if p>0
    local lap_mode  = n+1
    local lap_cmode = n+2
    Ord1 = 1
    h_lapl   = GLLapNLTerm(i,Ord,m,Vext,Ord1,lap_mode,lap_cmode,Lapg,LapCg)
    h_asymp .= h_asymp .+ h_lapl 
  end

  Ord1      = 1
  Ord2      = 1
  Ord3      = 1
  h_NL      = GLStdAsympNLTerm(i,Ord,m,Vext,Vext,Vext,Ord1,Ord2,Ord3)
  h_asymp  .= h_asymp .+ h_NL
  println("|NL|: $(norm(h_NL))")

  h_offd    = zeros(ComplexF64,N)
  for j in 1:(i-1)
    h_offd .= h_offd .+ SysMat_O2[i,j]*Y_O2[:,j]
  end
  h_asymp  .= h_asymp .- h_offd
  println("|h_offd|: $(norm(h_offd))")

  ω         = 0.0im
  for j in 1:Ord
    ii      = ind[j]+1
    ω       = ω + Khat[ii,ii] 
  end 

  DQ  = Matrix{ComplexF64}(I,Nby2,Nby2)
  CQ  = Matrix{ComplexF64}(I,Nby2,Nby2)
  for j in 1:n
    if abs(Khat[j,j] - ω) < 1.0e-12
      # println("Resonant ω: $(ω)")
      β           = W[:,j]'*(Bg2.*h_asymp)
      G_O2[j,i]   = β
      
      h_asymp    .= h_asymp .- β*V[:,j]
      DQ         .= DQ - V[ind1,j]*(W[ind1,j].*Bg)'
      CQ         .= CQ - V[ind2,j]*(W[ind2,j].*Bg)'
    end
  end  
  @views SEM1D.SEM_SetBC!(h_asymp[ind1],Inp.lbc,Inp.rbc)
  @views SEM1D.SEM_SetBC!(h_asymp[ind2],Inp.lbc,Inp.rbc)

  println("|r|: $(norm(h_asymp))")
  Y_O2[:,i]       = copy(h_asymp)

  Res1            = DQ*(ω*I - OPg)*DQ
  y_asymp         = view(Y_O2,ind1,i)
  r               = view(h_asymp,ind1)
  gmres!(y_asymp,Res1,r)

  Res2            = CQ*(ω*I - OPCg)*CQ
  y_asymp         = view(Y_O2,ind2,i)
  r               = view(h_asymp,ind2)
  gmres!(y_asymp,Res2,r)

end

# SysMat = BuildAsympSystem(Ord,m,Khat)

# # Third Order Asymptotic terms
# Ord   = 3
# Nt    = CenterManifold.NInteractionTerms(Ord,m)
# H3    = zeros(ComplexF64,N,Nt)
# 
# for i in 1:Nt
#   ind = CenterManifold.GetPolynomialIndices(i,Ord,m)
#   local i1  = ind[1]+1
#   local i2  = ind[2]+1
#   local i3  = ind[3]+1
# 
#   # println("($i1,$i2,$i3)")
# 
#   Ord1 = 1
#   Ord2 = 1
#   Ord3 = 1
#   h_asymp   = GLStdAsympNLTerm(i,Ord,m,Vext,Vext,Vext,Ord1,Ord2,Ord3)
#   H3[:,i]   = copy(h_asymp)
#   hnorm     = h_asymp'*(Bg2.*h_asymp)
#   # println("Index: $ind, Norm: $hnorm")
# end

println("Asymptotic System Done.")
















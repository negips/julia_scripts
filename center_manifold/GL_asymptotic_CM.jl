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

  Nt  = CenterManifold.NInteractionTerms(Ord,Nc)

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
            ind3 = [ind2[:]; (j-1)]
            
            k    = CenterManifold.GetIndexNumber(ind3,Nc)

            # println("Ind: $ind, k:$k")

            SysM[k,i] = SysM[k,i] + Khat[α,j]
          end     # j
        end       # if ...
      end         # i
    end           # ii
  end             # α

  return SysM
end  
#----------------------------------------------------------------------

# Second Order Asymptotic terms
Ord   = 2
Nt    = CenterManifold.NInteractionTerms(Ord,m)
H2    = zeros(ComplexF64,N,Nt)

for i in 1:Nt
  ind = CenterManifold.GetPolynomialIndices(i,Ord,m)
  local i1  = ind[1]+1
  local i2  = ind[2]+1

  # println("($i1,$i2)")

  local lap_mode  = n+1
  local lap_cmode = n+2
  Ord1 = 1
  h_asymp   = GLLapNLTerm(i,Ord,m,Vext,Ord1,lap_mode,lap_cmode,Lapg,LapCg)
  # @views SEM1D.SEM_SetBC!(h_asymp[ind1],Inp.lbc,Inp.rbc)
  # @views SEM1D.SEM_SetBC!(h_asymp[ind2],Inp.lbc,Inp.rbc)

  H2[:,i]   = copy(h_asymp)
  hnorm     = h_asymp'*(Bg2.*h_asymp)
  # println("Index: $ind, Norm: $hnorm")
  # if abs(hnorm)>1.0e-12
  #   leg = L"\mathfrak{R}(lapl"*"$i"*L")"
  #   ax2.plot(xg,real.(H2[ind1,i]),linewidth=1,linestyle="-", color=cm(m+i),label=leg)
  #   leg = L"\mathfrak{R}(lapl"*"$i"*L")"
  #   ax2.plot(xg,real.(H2[ind2,i]),linewidth=3,linestyle="--",color=cm(m+i),label=leg)
  # end  
end

SysMat = BuildAsympSystem(Ord,m,Khat)

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
















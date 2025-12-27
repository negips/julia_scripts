# Calculate the Center-Manifold Asymptotic vectors
include("../Module_CenterManifold/CenterManifold.jl")
#---------------------------------------------------------------------- 
function GLStdNLTerm(n0::Int,Ord0::Int,Nc::Int,MV1::AbstractMatrix{T},MV2::AbstractMatrix{T},MV3::AbstractMatrix{T},Ord1::Int,Ord2::Int,Ord3::Int,δ5) where {T}

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
function GLStdAsympNLTerm(n0::Int,Ord0::Int,Nc::Int,MV1::AbstractMatrix{T},MV2::AbstractMatrix{T},MV3::AbstractMatrix{T},Ord1::Int,Ord2::Int,Ord3::Int,δ5) where {T}

  if (Ord1 + Ord2 + Ord3 != Ord0) || Ord0 < 3
     h  = zeros(T,N)
     return h
  end

  if (Ord1 == Ord2) && (Ord1 == Ord3)

    h  =      GLStdNLTerm(n0,Ord0,Nc,MV1,MV2,MV3,Ord1,Ord2,Ord3,δ5)

  elseif (Ord1 == Ord2) && (Ord1 != Ord3)    

    h  =      GLStdNLTerm(n0,Ord0,Nc,MV1,MV2,MV3,Ord1,Ord2,Ord3,δ5)
    h .= h .+ GLStdNLTerm(n0,Ord0,Nc,MV1,MV3,MV2,Ord1,Ord3,Ord2,δ5)
    h .= h .+ GLStdNLTerm(n0,Ord0,Nc,MV3,MV1,MV2,Ord3,Ord1,Ord2,δ5)

  elseif (Ord1 == Ord3) && (Ord1 != Ord2)    

    h  =      GLStdNLTerm(n0,Ord0,Nc,MV1,MV2,MV3,Ord1,Ord2,Ord3,δ5)
    h .= h .+ GLStdNLTerm(n0,Ord0,Nc,MV2,MV1,MV3,Ord2,Ord1,Ord3,δ5)
    h .= h .+ GLStdNLTerm(n0,Ord0,Nc,MV3,MV1,MV2,Ord3,Ord1,Ord2,δ5)

  elseif (Ord2 == Ord3) && (Ord1 != Ord2)    

    h  =      GLStdNLTerm(n0,Ord0,Nc,MV1,MV2,MV3,Ord1,Ord2,Ord3,δ5)
    h .= h .+ GLStdNLTerm(n0,Ord0,Nc,MV2,MV1,MV3,Ord2,Ord1,Ord3,δ5)
    h .= h .+ GLStdNLTerm(n0,Ord0,Nc,MV2,MV3,MV1,Ord2,Ord3,Ord1,δ5)

  elseif (Ord1 != Ord2) && (Ord2 != Ord3) && (Ord3 != Ord1)

    h  =      GLStdNLTerm(n0,Ord0,Nc,MV1,MV2,MV3,Ord1,Ord2,Ord3,δ5)
    h .= h .+ GLStdNLTerm(n0,Ord0,Nc,MV1,MV3,MV2,Ord1,Ord3,Ord2,δ5)
    h .= h .+ GLStdNLTerm(n0,Ord0,Nc,MV2,MV1,MV3,Ord2,Ord1,Ord3,δ5)
    h .= h .+ GLStdNLTerm(n0,Ord0,Nc,MV2,MV3,MV1,Ord2,Ord3,Ord1,δ5)
    h .= h .+ GLStdNLTerm(n0,Ord0,Nc,MV3,MV1,MV2,Ord3,Ord1,Ord2,δ5)
    h .= h .+ GLStdNLTerm(n0,Ord0,Nc,MV3,MV2,MV1,Ord3,Ord2,Ord1,δ5)

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
function GLModePertTerm(n0::Int,Ord0::Int,MV1::AbstractMatrix{T1},Ord1::Int,V::AbstractMatrix{T1},W::AbstractMatrix{T1},B::AbstractVector{T2},PertModes::AbstractVector{Int}) where {T1,T2<:Number}

  ngs       = 2
  N,Nt1     = size(MV1)                # Vector Length, No of Modes
  Nc        = size(V,2)                # No. of Critical Modes
  npmodes   = length(PertModes) 

  h         = zeros(T1,N)
  vtmp      = zeros(T1,N)
  if (Ord1 + 1 != Ord0) 
    return h
  end
  if Ord0<2
    return h
  end
  if npmodes == 0
    return h
  end  

  # Index we are trying to match
  ind0      = CenterManifold.GetPolynomialIndices(n0,Ord0,Nc) .+ 1

  for i in 1:Nt1
    ind1        = CenterManifold.GetPolynomialIndices(i,Ord1,Nc) .+ 1
    for j in 1:npmodes
      # Mode Perturbation nonlinearity σjVj⋅(B⋅Wj)'
      pm = PertModes[j]
      if pm>0
        ind_total   = [ind1[:]; pm]
        sort!(ind_total)
        if (ind_total == ind0)    # Matching polynomials
          di  = 1
          si  = (i-1)*N + 1
          copyto!(vtmp,di,MV1,si,N)
          β   = T1(0)
          for k in 1:ngs
            α       = W[:,j]'*(B.*vtmp)
            β       = β + α
            vtmp   .= vtmp .- α*V[:,j]
          end
          h .= h .+ β*V[:,j]
        end       # ind_total == ind0
      end         # pm>0
    end           # j in 1:nmodes
  end             # i in 1:Nt1

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
            ind3        = CenterManifold.GetPolynomialIndices(j,OrdG,Nc)
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

















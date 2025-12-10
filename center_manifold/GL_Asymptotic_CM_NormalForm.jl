# Testing the module

# include("Module_SEM1D/SEM1D.jl")
# #using .SEM1D
# 
# include("Module_StepperArnoldi/StepperArnoldi.jl")
# #using .StepperArnoldi
# 
include("Module_CenterManifold/CenterManifold.jl")
# #using .CenterManifold

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

  # println(ind0)
  # println(norm(h))

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
#---------------------------------------------------------------------- 

h3      = figure(num=3,figsize=Grh.figsz2)
ax3     = gca()
cm2     = get_cmap("tab10")
plno    = -1
verbose = false
plforc  = true


# Second Order Asymptotic terms
#-------------------------------------------------- 
Ord         = 2
Nt          = CenterManifold.NInteractionTerms(Ord,m)
Y_O2        = zeros(ComplexF64,N,Nt)
SysMat_O2   = BuildAsympSystem(Ord,m,Khat)

# Reduced Matrix terms
G_O2        = zeros(ComplexF64,m,Nt)
println("Solving for Ord=$(Ord)")
for i in 1:Nt

  global plno

  ind       = CenterManifold.GetPolynomialIndices(i,Ord,m)
  ω         = 0.0im
  for j in 1:Ord
    ii      = ind[j]+1
    ω       = ω + Khat[ii,ii] 
  end 

  if (verbose)
    println("Solving for $(ind .+ 1), ω=$(ω)")
  end  

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
  δ5        = δ[5]
  h_NL      = GLStdAsympNLTerm(i,Ord,m,Vext,Vext,Vext,Ord1,Ord2,Ord3,δ5)
  h_asymp  .= h_asymp .+ h_NL
  # println("|NL|: $(norm(h_NL))")

  h_offd    = zeros(ComplexF64,N)
  for j in 1:(i-1)
    h_offd .= h_offd .+ SysMat_O2[i,j]*Y_O2[:,j]
  end
  h_asymp  .= h_asymp .- h_offd
  # println("|h_offd|: $(norm(h_offd))")

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

  # println("|r|: $(norm(h_asymp))")
  Y_O2[:,i]       = copy(h_asymp)

  Res1            = DQ*(ω*I - OPg)*DQ
  y_asymp         = view(Y_O2,ind1,i)
  r               = view(h_asymp,ind1)
  gmres!(y_asymp,Res1,r)
  y1norm = sqrt(abs(y_asymp'*(Bg.*y_asymp))) 
  # plotting
  if y1norm>1.0e-12
    plno    = plno + 1
    ijk = "_{$(ind[1]+1)"
    for k in 2:Ord
      ijk = ijk*",$(ind[k]+1)"
    end
    ijk = ijk*"}"
    # leg = L"\mathfrak{R}(y%$ijk)"
    # ax3.plot(xg,real.(y_asymp),linewidth=2,linestyle="-", color=cm2(plno),label=leg)
    # leg = L"\mathfrak{Im}(y%$ijk)"
    # ax3.plot(xg,imag.(y_asymp),linewidth=2,linestyle="--",color=cm2(plno),label=leg)

    leg = L"|y%$ijk|"
    ax3.plot(xg,abs.(y_asymp),linewidth=2,linestyle="-", color=cm2(plno),label=leg)

    # leg = L"|h%$ijk|"
    # ax3.plot(xg,abs.(h_asymp[ind1]),linewidth=1,linestyle="--", color=cm2(plno),label=leg)

    # leg = L"|hNL%$ijk|"
    # ax3.plot(xg,abs.(h_NL[ind1]),linewidth=1,linestyle="-.", color=cm2(plno),label=leg)
    # println("|NL1|: $(norm(h_NL[ind1]))")

    # leg = L"|hLap%$ijk|"
    # ax3.plot(xg,abs.(h_lapl[ind1]),linewidth=1,linestyle="-.", color=cm2(plno),label=leg)
  
  end  

  Res2            = CQ*(ω*I - OPCg)*CQ
  y_asymp         = view(Y_O2,ind2,i)
  r               = view(h_asymp,ind2)
  gmres!(y_asymp,Res2,r)
  y2norm = sqrt(abs(y_asymp'*(Bg.*y_asymp))) 
  # plotting no
  if y2norm>1.0e-12
    plno    = plno + 1
    ijk = "_{$(ind[1]+1)"
    for k in 2:Ord
      ijk = ijk*",$(ind[k]+1)"
    end
    ijk = ijk*"}"
    # leg = L"\mathfrak{R}(y%$ijk*)"
    # ax3.plot(xg,real.(y_asymp),linewidth=2,linestyle="-", color=cm2(plno),label=leg)
    # leg = L"\mathfrak{Im}(y%$ijk*)"
    # ax3.plot(xg,imag.(y_asymp),linewidth=2,linestyle="--",color=cm2(plno),label=leg)

    # leg = L"|y%$ijk*|"
    # ax3.plot(xg,abs.(y_asymp),linewidth=2,linestyle="--",color=cm2(plno),label=leg)
  end  

end
ax3.legend(ncol=3,loc="upper right",fontsize=Grh.lgfs)

#-------------------------------------------------- 


h4    = figure(num=4,figsize=Grh.figsz2)
ax4   = gca()
plno  = -1

# Third Order Asymptotic terms
#-------------------------------------------------- 
Ord         = 3
Nt          = CenterManifold.NInteractionTerms(Ord,m)
Y_O3        = zeros(ComplexF64,N,Nt)
SysMat_O3   = BuildAsympSystem(Ord,m,Khat)

# Reduced Matrix terms
G_O3        = zeros(ComplexF64,m,Nt)
println("Solving for Ord=$(Ord)")
for i in 1:Nt

  global plno

  ind       = CenterManifold.GetPolynomialIndices(i,Ord,m)
  ω         = 0.0im
  for j in 1:Ord
    ii      = ind[j]+1
    ω       = ω + Khat[ii,ii] 
  end 

  if (verbose)
    println("Solving for $(ind .+ 1), ω=$(ω)")
  end  

  h_asymp   = zeros(ComplexF64,N)
  if p>0
    local lap_mode  = n+1
    local lap_cmode = n+2
    OrdY     = 2
    h_lapl   = GLLapNLTerm(i,Ord,m,Y_O2,OrdY,lap_mode,lap_cmode,Lapg,LapCg)
    h_asymp .= h_asymp .+ h_lapl 
  end

  Ord1      = 1
  Ord2      = 1
  Ord3      = 1
  δ5        = δ[5]
  h_NL      = GLStdAsympNLTerm(i,Ord,m,Vext,Vext,Vext,Ord1,Ord2,Ord3,δ5)
  h_asymp  .= h_asymp .+ h_NL

  h_offd    = zeros(ComplexF64,N)
  for j in 1:(i-1)
    h_offd .= h_offd .+ SysMat_O3[i,j]*Y_O3[:,j]
  end
  h_asymp  .= h_asymp .- h_offd

  OrdY      = 2
  OrdG      = 2
  h_low     = BuildLowOrder(i,Ord,m,Y_O2,G_O2,OrdY,OrdG)
  h_asymp  .= h_asymp .- h_low
  #println("|r_low|: $(norm(h_low))")

  DQ  = Matrix{ComplexF64}(I,Nby2,Nby2)
  CQ  = Matrix{ComplexF64}(I,Nby2,Nby2)
  for j in 1:n
    if abs(Khat[j,j] - ω) < 1.0e-12
      # println("Resonant ω: $(ω)")
      β           = W[:,j]'*(Bg2.*h_asymp)
      G_O3[j,i]   = β
      
      h_asymp    .= h_asymp .- β*V[:,j]
      DQ         .= DQ - V[ind1,j]*(W[ind1,j].*Bg)'
      CQ         .= CQ - V[ind2,j]*(W[ind2,j].*Bg)'
    end
  end  
  @views SEM1D.SEM_SetBC!(h_asymp[ind1],Inp.lbc,Inp.rbc)
  @views SEM1D.SEM_SetBC!(h_asymp[ind2],Inp.lbc,Inp.rbc)

  Y_O3[:,i]       = copy(h_asymp)

  Res1            = DQ*(ω*I - OPg)*DQ
  y_asymp         = view(Y_O3,ind1,i)
  r               = view(h_asymp,ind1)
  gmres!(y_asymp,Res1,r)
  y1norm = sqrt(abs(y_asymp'*(Bg.*y_asymp))) 
  # plotting
  if y1norm>1.0e-12
    plno    = plno + 1
    ijk = "_{$(ind[1]+1)"
    for k in 2:Ord
      ijk = ijk*",$(ind[k]+1)"
    end
    ijk = ijk*"}"
    # leg = L"\mathfrak{R}(y%$ijk)"
    # ax4.plot(xg,real.(y_asymp),linewidth=2,linestyle="-", color=cm2(plno),label=leg)
    # leg = L"\mathfrak{Im}(y%$ijk)"
    # ax4.plot(xg,imag.(y_asymp),linewidth=2,linestyle="--",color=cm2(plno),label=leg)
    
    leg = L"|y%$ijk|"
    ax4.plot(xg,abs.(y_asymp),linewidth=2,linestyle="-", color=cm2(plno),label=leg)
  end  

  Res2            = CQ*(ω*I - OPCg)*CQ
  y_asymp         = view(Y_O3,ind2,i)
  r               = view(h_asymp,ind2)
  gmres!(y_asymp,Res2,r)
  y2norm = sqrt(abs(y_asymp'*(Bg.*y_asymp))) 
  # plotting no
  if y2norm>1.0e-12
    plno    = plno + 1
    ijk = "_{$(ind[1]+1)"
    for k in 2:Ord
      ijk = ijk*",$(ind[k]+1)"
    end
    ijk = ijk*"}"
    # leg = L"\mathfrak{R}(y%$ijk*)"
    # ax4.plot(xg,real.(y_asymp),linewidth=2,linestyle="-", color=cm2(plno),label=leg)
    # leg = L"\mathfrak{Im}(y%$ijk*)"
    # ax4.plot(xg,imag.(y_asymp),linewidth=2,linestyle="--",color=cm2(plno),label=leg)

    # leg = L"|y%$ijk*|"
    # ax4.plot(xg,abs.(y_asymp),linewidth=2,linestyle="--",color=cm2(plno),label=leg)
  end  

end

ax4.legend(ncol=3,loc="upper right",fontsize=Grh.lgfs)

println("Asymptotic System Done.")
















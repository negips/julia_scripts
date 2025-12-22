# Perturbation of the Ginzburg Landau Eigenmodes
#---------------------------------------------------------------------- 
function GLSelectPertModes(DOut,AOut,modeselect::Vector{Int})

      T1          = eltype(DOut.evecs)
      T2          = eltype(DOut.evals)
      N           = size(DOut.evecs,1)
      N2          = N*2
      nmode       = length(modeselect)
      npert       = nmode*2
      Vpert       = zeros(T1,N2,npert)
      Wpert       = zeros(T1,N2,npert)
      ind1        = 1:N
      ind2        = N+1:N2
      λpert       = zeros(T2,npert) 
      λpertA      = zeros(T2,npert)     # Adjoint

      λ           = copy(DOut.evals)
      sind        = sortperm(real.(λ),rev=true)       # Sorted (real part) indices
      si          = sind[modeselect]                  # Selected indices
      for i in 1:nmode
        j1              = (i-1)*2 + 1
        j2              = j1 + 1
        ii              = si[i]
        Vpert[ind1,j1]  = copy(DOut.evecs[:,ii])
        Vpert[ind2,j2]  = conj.(Vpert[ind1,j1])
        λpert[j1]       = DOut.evals[ii]
        λpert[j2]       = conj(λpert[j1])
      end

      λ           = copy(AOut.evals)
      sind        = sortperm(real.(λ),rev=true)    # Sorted (real part) indices
      si          = sind[modeselect]               # Selected indices
      for i in 1:nmode
        j1              = (i-1)*2 + 1
        j2              = j1 + 1
        ii              = si[i]
        Wpert[ind1,j1]  = copy(AOut.evecs[:,ii])
        Wpert[ind2,j2]  = conj.(Wpert[ind1,j1])
        λpertA[j1]      = AOut.evals[ii]
        λpertA[j2]      = conj(λpertA[j1])
      end

      return λpert,λpertA,Vpert,Wpert
end  
#---------------------------------------------------------------------- 
function GLModePertTerm(n0::Int,Ord0::Int,Nc::Int,MV1::AbstractMatrix{T},Ord1::Int,V::AbstractMatrix{T},W::AbstractMatrix{T},B::AbstractVector{T1},PertModes::AbstractVector{Int}) where {T,T1<:Number}

  ngs  = 2
  ind0 = CenterManifold.GetPolynomialIndices(n0,Ord0,Nc) .+ 1

  Nt1  = CenterManifold.NInteractionTerms(Ord1,Nc)

  N,c  = size(MV1)
  
  nmodes = length(PertModes) 

  h    = zeros(T,N)
  vtmp = zeros(T,N)
  if (Ord1 + 1 != Ord0) 
    return h
  end
  if Ord0<2
    return h
  end
  if nmodes == 0
    return h
  end  


  for i in 1:Nt1
    for j in PertModes
      # Mode Perturbation nonlinearity σjVj⋅(B⋅Wj)'
      ind1        = CenterManifold.GetPolynomialIndices(i,Ord1,Nc) .+ 1
      ind_total   = [ind1[:]; PertModes[j]]
      sort!(ind_total)
      if (ind_total == ind0)    # Matching polynomials
        di  = 1
        si  = (i-1)*N + 1
        copyto!(vtmp,di,MV1,si,N)
        β   = T(0)
        for k in 1:ngs
          α       = Wpert[:,j]*(B.*vtmp)
          β       = β + α
          vtmp   .= vtmp .- α*Vpert[:,j]
        end
        h .= h .+ β*V[:,j]
      end         # ind_total == ind0
    end           # j in 1:nmodes
  end             # i in 1:Nt1

  return h
end
#---------------------------------------------------------------------- 

ifmodepert = false

if (ifmodepert)
  modeselect                    = [2]
  nmode                         = length(modeselect)
  λpert,λpertA,Vpert,Wpert      = GLSelectPertModes(ArnDir,ArnAdj,modeselect)
  npert                         = length(λpert)
  λnew                          = zeros(eltype(λpert),npert)
  dλ                            = zeros(eltype(λpert),npert)
  
  # Renormalize
  j0    = renormalization_index(xg)
  for i in 1:nmode
    j1 = (i-1)*2 + 1
    j2 = j1 + 1
    @views renormalize_evec!(Vpert[ind1,j1],j0)
    @views renormalize_evecs!(Vpert[ind1,j1],Wpert[ind1,j1],Bg)
 
    Vpert[ind2,j2] = conj.(Vpert[ind1,j1])
    @views renormalize_evecs!(Vpert[ind2,j2],Wpert[ind2,j2],Bg)
  end  
  
  for i in 1:length(λpert)
    λnew[i] = imag(λpert[i])*im       # Map to real axis
    dλ[i]   = (λnew[i] - λpert[i])*1
  end  

  nsys        = n + npert
  VSys        = [V Vpert]
  WSys        = [W Wpert]
  λSys        = [λc; λnew]
  σ           = [zeros(ComplexF64,n);dλ]
  PertModes   = zeros(Int,nsys)
  PertModes[n+1:nsys] = nsys+1:nsys+npert

  VSys1       = VSys[ind1,:]
  VSys2       = VSys[ind2,:]
  WSys1       = WSys[ind1,:]
  WSys2       = WSys[ind2,:]

  # V1          = VSys[ind1,1:2:nsys]
  # V2          = VSys[ind2,2:2:nsys]
  # W1          = WSys[ind1,1:2:nsys]
  # W2          = WSys[ind2,2:2:nsys]
  # σ1          = σ[1:2:nsys]
  # σ2          = σ[2:2:nsys]
  # PModes1     = PertModes[1:2:nsys]
  # PModes2     = PertModes[2:2:nsys]

  # V1          = Vpert[ind1,1:2:npert]
  # V2          = Vpert[ind2,2:2:npert]
  # W1          = Wpert[ind1,1:2:npert]
  # W2          = Wpert[ind2,2:2:npert]
  # σ1          = dλ[1:2:npert]
  # σ2          = dλ[2:2:npert]
 
  # New Linear Operator
  # NewOPg(x)   = PertLv(x,OPg,Bg,VSys1,WSys1,σ)      
  # NewOPCg(x)  = PertLv(x,OPCg,Bg,VSys2,WSys2,σ)

  # NewEOPg(x,y,z)  = ExtPertLv(x,OPg,Bg,VSys1,WSys1,σ,y,z)      
  # NewEOPCg(x,y,z) = ExtPertLv(x,OPCg,Bg,VSys2,WSys2,σ,y,z)

  # ArnDirNew   = StepperArnoldi.StepArnOP(NewOPg,Bg,StpInp,ArnInp,Inp.lbc,Inp.rbc)

else

  nmode       = 0
  npert       = 0
  nsys        = n + npert
  λSys        = copy(λc)
  σ           = zeros(ComplexF64,nsys)

  PertModes   = zeros(Int,nsys)

  VSys        = copy(V)
  WSys        = copy(W)

  VSys1       = VSys[ind1,:]
  VSys2       = VSys[ind2,:]
  WSys1       = WSys[ind1,:]
  WSys2       = WSys[ind2,:]

  # NewOPg      = OPg
  # NewOPCg     = OPCg

end  

println("Eigenmode Perturbation Done.")
















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
    λ1      = sign(imag(λpert[i]))*im       # Map to real axis
    dλ[i]   = (λ1 - λpert[i])*1
    λnew[i] = λpert[i] + dλ[i] 
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

  for i in 1:nmode
    j1      = (i-1)*2 + 1
    j2      = j1 + 1
    vtmp    = Vpert[ind1,j1]
    wtmp    = Wpert[ind1,j1]
 
    j       = n + j1

    ax2.plot(xg,real.(vtmp),linewidth=2,linestyle="-", color=cm(j-1),label=L"\mathfrak{R}(ϕ_{%$j})")
    ax2.plot(xg,imag.(vtmp),linewidth=2,linestyle="--",color=cm(j-1),label=L"\mathfrak{Im}(ϕ_{%$j})")
  
    ax2.plot(xg,real.(wtmp),linewidth=1,linestyle="-", color=cm(j-1+10),label=L"\mathfrak{R}(χ_{%$j})")
    ax2.plot(xg,imag.(wtmp),linewidth=1,linestyle="--",color=cm(j-1+10),label=L"\mathfrak{Im}(χ_{%$j})")
  end
  ax2.legend(ncols=2,fontsize=Grh.lgfs)

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

end  

println("Eigenmode Perturbation Done.")
















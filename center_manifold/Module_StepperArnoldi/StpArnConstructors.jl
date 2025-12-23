function ExtendedModes(EM::Vector{ExtendedMode})

  T  = typeof(EM[1].λe)       # Type
  ne = length(EM)             # No of Extended Modes
  N  = length(EM[1].ve)       # Original system Mode length
  nc = size(EM[1].z,2)        # Critical space

  Γ  = zeros(T,nc,ne)
  Fe = zeros(T,N,ne)
  Ve = zeros(T,N,ne)
  We = zeros(T,N,ne)
  Z  = zeros(T,ne,nc)
  λe = zeros(T,ne)

  for i in 1:ne
    λe[i] = EM[i].λe
    for j in 1:nc
      Γ[j,i] = EM[i].Γ[j]
      Z[i,j] = EM[i].z[j]
    end
    di = (i-1)*N+1
    si = 1
    copyto!(Fe,di,EM[i].fe,si,N)
    copyto!(Ve,di,EM[i].ve,si,N)
    copyto!(We,di,EM[i].we,si,N)
  end  

  ExModes = ExtendedModes(λe,Fe,Γ,Ve,We,Z)

  return ExModes
end
#----------------------------------------------------------------------       


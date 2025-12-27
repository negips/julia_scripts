# Calculate the Center-Manifold Asymptotic vectors
include("GLAsymptoticCM.jl")
#---------------------------------------------------------------------- 
function plotAsymptoticVectors(ax,xg::AbstractVector{T1},Y::AbstractMatrix{T2},Ord::Int,Nc::Int,lgfs) where {T1,T2<:Number}

  nt   = CenterManifold.NInteractionTerms(Ord,Nc)
  plno = -1
  cmap = get_cmap("tab20")

  for i in 1:nt
    vnorm = norm(Y[:,i])
    # Plot Mode
    if vnorm > 0.0
   
      ind = CenterManifold.GetPolynomialIndices(i,Ord,Nc) .+ 1
      ijk = "_{$(ind[1])"
      for k in 2:Ord
        ijk = ijk*",$(ind[k])"
      end
      ijk = ijk*"}"
     
      j = i
      plno = plno + 1 
      ax.plot(xg,real.(Y[:,i]),linewidth=2,linestyle="-", color=cmap(plno),label=L"\mathfrak{R}(y_{%$ijk})")
      ax.plot(xg,imag.(Y[:,i]),linewidth=2,linestyle="--",color=cmap(plno),label=L"\mathfrak{Im}(y_{%$ijk})")
    end
  end
  plno = plno + 1

  if plno<=4
    ax.legend(fontsize=lgfs,loc="upper right")
  elseif plno>4 && plno<=12
    ax.legend(ncols=3,fontsize=lgfs,loc="upper right")
  elseif plno>12 && plno<=20
    ax.legend(ncols=4,fontsize=lgfs,loc="upper right")
  elseif plno>20 && plno<=30
    ax.legend(ncols=5,fontsize=lgfs,loc="upper right")
  end

  return nothing
end  
#---------------------------------------------------------------------- 
resmodeplot = true

MaxOrd      = 5

σext  = zeros(ComplexF64,m)
copyto!(σext,1,σ,1,nsys)
PertModesExt      = zeros(Int64,m)
copyto!(PertModesExt,1,PertModes,1,nsys)

# Ord = 2
SM2,G2,Y2,F2,CMmodes2 = GL2AsymptoticCM_Ord2(OPg,OPCg,Bg,δ,Khat,Vext,Wext,σext,PertModesExt,nsys,p,h,Inp.lbc,Inp.rbc,BiLapg,BiLapCg,ifnormal);
if (resmodeplot)
  h4    = figure(num=4,figsize=Grh.figsz2)
  ax4   = gca()
  ax4.cla()
  @views plotAsymptoticVectors(ax4,xg,Y2[ind1,:],2,m,Grh.lgfs)
end

# Ord = 3
SM3,G3,Y3,F3,CMmodes3 = GL2AsymptoticCM_Ord3(OPg,OPCg,Bg,δ,Khat,Vext,Wext,σext,PertModesExt,nsys,p,h,G2,Y2,Inp.lbc,Inp.rbc,BiLapg,BiLapCg,ifnormal);
if (resmodeplot)
  h5    = figure(num=5,figsize=Grh.figsz2)
  ax5   = gca()
  ax5.cla()
  @views plotAsymptoticVectors(ax5,xg,Y3[ind1,:],3,m,Grh.lgfs)
end

if MaxOrd>=4
  # Ord = 4
  SM4,G4,Y4,F4,CMmodes4 = GL2AsymptoticCM_Ord4(OPg,OPCg,Bg,δ,Khat,Vext,Wext,σext,PertModesExt,nsys,p,h,G2,Y2,G3,Y3,Inp.lbc,Inp.rbc,BiLapg,BiLapCg,ifnormal);
  if (resmodeplot)
    h6    = figure(num=6,figsize=Grh.figsz2)
    ax6   = gca()
    ax6.cla()
    @views plotAsymptoticVectors(ax6,xg,Y4[ind1,:],4,m,Grh.lgfs)
  end
end  

if MaxOrd>=5
  # Ord = 5
  SM5,G5,Y5,F5,CMmodes5 = GL2AsymptoticCM_Ord5(OPg,OPCg,Bg,δ,Khat,Vext,Wext,σext,PertModesExt,nsys,p,h,G2,Y2,G3,Y3,G4,Y4,Inp.lbc,Inp.rbc,BiLapg,BiLapCg,ifnormal);
  if (resmodeplot)
    h7    = figure(num=7,figsize=Grh.figsz2)
    ax7   = gca()
    ax7.cla()
    @views plotAsymptoticVectors(ax7,xg,Y5[ind1,:],5,m,Grh.lgfs)
  end
end  

println("Asymptotic System (Normal Form) Done to Ord=$MaxOrd.")
















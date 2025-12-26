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

σext  = zeros(ComplexF64,m)
copyto!(σext,1,σ,1,nsys)
PertModesExt      = zeros(Int64,m)
copyto!(PertModesExt,1,PertModes,1,nsys)


SM2,G2,Y2,F2,CMmodes2 = GL2AsymptoticCM_Ord2(OPg,OPCg,Bg,δ,Khat,Vext,Wext,σext,PertModesExt,nsys,p,h,Inp.lbc,Inp.rbc,BiLapg,BiLapCg,ifnormal);
if (resmodeplot)
  h4    = figure(num=4,figsize=Grh.figsz2)
  ax4   = gca()
  ax4.cla()
  
  @views plotAsymptoticVectors(ax4,xg,Y2[ind1,:],2,m,Grh.lgfs)
end

SM3,G3,Y3,F3,CMmodes3 = GL2AsymptoticCM_Ord3(OPg,OPCg,Bg,δ,Khat,Vext,Wext,σext,PertModesExt,nsys,p,h,G2,Y2,Inp.lbc,Inp.rbc,BiLapg,BiLapCg,ifnormal);
 
if (resmodeplot)
  h5    = figure(num=5,figsize=Grh.figsz2)
  ax5   = gca()
  ax5.cla()

  @views plotAsymptoticVectors(ax5,xg,Y3[ind1,:],3,m,Grh.lgfs)
end

println("Asymptotic System (Normal Form) Done.")
















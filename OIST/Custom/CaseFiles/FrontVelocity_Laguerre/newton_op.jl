function newton_Lop(x,C0,gradf,Lap,Conv,B)

  Ax      = (Lap*x .+ C0*Conv*x .+ B.*gradf.*x);
  Ax[1]   = 0.0
  Ax[end] = 0.0

  # v         = Vector(Conv[end,:])
  # Ax[end-1] = v'*x   # Putting ∇x == 0 BC at second last point

  return Ax
end 
#---------------------------------------------------------------------- 
function newton_resid(x,C0,fop,Lap,Conv,B,bcl,bcr)

  resid      = -(Lap*x .+ C0*Conv*x .+ B.*fop(x));
  resid[1]   = bcl
  resid[end] = bcr

  # v         = Vector(Conv[end,:])
  # resid[end-1] = -v'*x   # Putting ∇x at second last point

  return resid
end
#----------------------------------------------------------------------
function newton_Lop_spectral(x,Basis,C0,gradf,Lap,Conv,M2N,N2M)

  Ax      = (Lap*x .+ C0*Conv*x .+ M2N'*diagm(Basis.orthoweights.*Basis.weights)*diagm(gradf)*(M2N*x) );

  # Operator Changed for Left Boundary condition
  Ax[1]   = (M2N[1,:])'*x
  Ax[end] = (M2N[end,:])'*x

  #Ax[1]   = 0.0
  #Ax[end] = 0.0

  # v         = Vector(Conv[end,:])
  # Ax[end-1] = v'*x   # Putting ∇x == 0 BC at second last point

  return Ax
end 
#---------------------------------------------------------------------- 
function newton_resid_spectral(x,Basis,C0,fop,Lap,Conv,M2N,N2M,bcl,bcr)

  resid      = -(Lap*x .+ C0*Conv*x .+ M2N'*diagm(Basis.orthoweights.*Basis.weights)*fop(M2N*x));
  resid[1]   = -(bcl - (M2N[1,:])'*x)
  resid[end] = -(bcr - (M2N[end,:])'*x)
 
  #resid[end] = bcr

  # v         = Vector(Conv[end,:])
  # resid[end-1] = -v'*x   # Putting ∇x at second last point

  return resid
end
#----------------------------------------------------------------------




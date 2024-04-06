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
function newton_resid_laguerre(x,C0,fop,Lap,Conv,B,bcl,bcr)

  resid      = -(Lap*x .+ C0*Conv*x .+ B.*fop(x));
  resid[1]   = bcl
  resid[end] = bcr

  # v         = Vector(Conv[end,:])
  # resid[end-1] = -v'*x   # Putting ∇x at second last point

  return resid
end
#----------------------------------------------------------------------




# Write out an hdf5 file

using HDF5

function write_params(group,f,g,h,pars,parsS,Da,ϵ,ν,ndigits)

  # f-parameters
  α1        = pars.fcy[3]
  α3        = pars.fcx[1]

  f00(y)    = f(0.0,y,0.0)
  roots_f   = find_zeros(f00,-10.0,10.0)

  j = 0
  α2 = -999.9
  for i in 1:length(roots_f)
    r = roots_f[i]
    if (r > 1.1)
      α2    = r
    end
  end

  # g-parameters
  β0        = pars.gc0
  β1        = pars.gcy[1]
  β2        = pars.gcx[1]
  R         = sqrt(β1^2 + β2^2)
  ϕ         = atan(β2,β1)

  # h-params
  γ0        = parsS.fc0
  γ1        = parsS.fcx[1]
  
  γ2        = parsS.fcy[3]
  abar      = 0.00
  h0(y)     = h(abar,y)
  hroots_0  = find_zeros(h0,-10.0,10.0)
  γ1        = γ1/abs(γ2)

  @printf "f-params\n"
  write_dataset(group,"α1",round(α1,digits=ndigits)) 
  write_dataset(group,"α2",round(α2,digits=ndigits)) 
  write_dataset(group,"α3",round(α3,digits=ndigits)) 
  write_dataset(group,"D",Da) 
  write_dataset(group,"ϵ",ϵ) 

  @printf "g-params\n"
  write_dataset(group,"β0",round(β0,digits=ndigits)) 
  write_dataset(group,"β1",round(β1,digits=ndigits)) 
  write_dataset(group,"β2",round(β2,digits=ndigits)) 

  @printf "h-params\n"
  write_dataset(group,"γ1",round(γ1,digits=ndigits)) 
  # write_dataset(group,"γ2",round(γ2,digits=ndigits)) 
  write_dataset(group,"ν",ν/abs(γ2))

  return nothing
end 
#---------------------------------------------------------------------- 

fname = "wedges.h5"
fid   = h5open(fname, "w")
# parameter data        
g1 = create_group(fid,"Params") 
#g2 = create_group(fid,"Data")        

# # Parameters
ndigits = 4
Da      = γa/ϵ
write_params(g1,F,G,λdot1,pars,parsS,Da,ϵ,ν,ndigits)

# write_dataset(g1,"ndim",re2.hdr.ldim)
# write_dataset(g1,"nelgv",re2.hdr.nelgv)
# write_dataset(g1,"nelgt",re2.hdr.nelgt)
# write_dataset(g1,"wdsize",re2.hdr.wdsize)
# write_dataset(g1,"ncurve",re2.ncurve)
# 
# # Data
# write_dataset(g2,"xc",re2.xc)
# write_dataset(g2,"yc",re2.yc)
# write_dataset(g2,"zc",re2.zc)
# write_dataset(g2,"bcs",re2.cbl)
# write_dataset(g2,"bcparams",re2.bl)
# if (re2.ncurve>0)
#   write_dataset(g2,"curveieg",re2.curveieg)
#   write_dataset(g2,"curveiside",re2.curveiside)
#   write_dataset(g2,"curveparam",re2.curveparam)
#   write_dataset(g2,"curvetype",re2.curvetype)
# end  
# 
# # .ma2 
# h  = create_group(fid,"Ma2")
# h1 = create_group(h,"Params")
# h2 = create_group(h,"Data")
# 
# # Parameters
# write_dataset(h1,"d2",map.hdr.d2)
# write_dataset(h1,"depth",map.hdr.depth)
# write_dataset(h1,"nactive",map.hdr.nactive)
# write_dataset(h1,"nel",map.hdr.nel)
# write_dataset(h1,"noutflow",map.hdr.noutflow)
# write_dataset(h1,"npts",map.hdr.npts)
# write_dataset(h1,"nrank",map.hdr.nrank)
# write_dataset(h1,"ifsorted",ifsorted)
# 
# # Data
# write_dataset(h2,"pmap",map.pmap)
# write_dataset(h2,"vmap",map.vmap)
close(fid)


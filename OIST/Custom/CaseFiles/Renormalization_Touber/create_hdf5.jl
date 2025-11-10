# Write out an hdf5 file

using HDF5

function write_params(group,f,g,h,pars,parsS,Da,β3,ϵ,ν,Asen,Aeq,ndigits,ifλ)

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
  if (ifλ)
    write_dataset(group,"β3",round(β3,digits=ndigits)) 
  end  

  if (ifλ)
    @printf "h-params\n"
    write_dataset(group,"γ1",round(γ1,digits=ndigits)) 
    # write_dataset(group,"γ2",round(γ2,digits=ndigits)) 
    write_dataset(group,"ν",ν/abs(γ2))
    write_dataset(group,"δ",Asen)
    write_dataset(group,"Aeq",Aeq)
  end  

  return nothing
end 
#---------------------------------------------------------------------- 
function write_data(group,x,time,B_hist,A_hist,λ_hist,ifλ)

  # x
  write_dataset(group,"x",x) 
  # Time
  write_dataset(group,"time",time) 
  # Inhibitor
  write_dataset(group,"B_hist",B_hist)
  # Activator
  write_dataset(group,"A_hist",A_hist)

  if (ifλ)
    write_dataset(group,"λ_hist",λ_hist)
  end  

  return nothing
end 
#---------------------------------------------------------------------- 

fid     = h5open(fnameh5, "w")
# parameter data        
g1 = create_group(fid,"Params") 
#g2 = create_group(fid,"Data")        

# # Parameters
ndigits = 4
Da      = γa/ϵ
write_params(g1,F,G,λdot1,pars,parsS,Da,β3,ϵ,ν,Asen,Aeq,ndigits,ifλ)

g2 = create_group(fid,"Data")
write_data(g2,Geom.xm1[:],Thist,fldhist[:,:,1],fldhist[:,:,2],γhist,ifλ)

@printf("%s HDF file saved\n",fnameh5)

close(fid)



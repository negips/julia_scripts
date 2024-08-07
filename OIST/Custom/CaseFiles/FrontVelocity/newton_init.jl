
println("Newton Initialization")

agauss      = 0.0*Geom.xm1[:]
if ngauss == 1
  agauss      = exp.(-((Geom.xm1[:] .- x0)/σg).^2)
else  
  for i in 1:ngauss
    global agauss
    agauss    = agauss .+ ampgauss[i]*exp.(-((Geom.xm1[:] .- x0gauss[i])/σg).^2)
  end
end  

#agauss      = exp.(-((Geom.xm1[:] .- x0)/σg).^2)# .*(sign.(Geom.xm1[:] .- x0))
k0          = 7
for ci in CartesianIndices(Geom.xm1)
  if isinf(Geom.xm1[ci])
    i,j = Tuple(ci)
    if i==1 && j==1
      Geom.xm1[ci] = -999.999
    else
      Geom.xm1[ci] =  999.999
    end
  end
end
asin        = sin.(2.0*π/(xe-xs)*k0*Geom.xm1[:])
acos        = cos.(2.0*π/(xe-xs)*k0*Geom.xm1[:])

ainit       = vimultg.*(QT*asin)
binit       = vimultg.*(QT*acos)

fld         = zeros(VT,ndof,nflds)

for j in 1:nflds
  fld[:,j]    = Amp0[j]*ainit .+ Off0[j] .+ σ0i[j]*(rand(ndof) .- 0.5)
end  


cm          = get_cmap("tab10");
rgba0       = cm(0); 
rgba1       = cm(1); 
rgba2       = cm(2);

h2          = figure(num=2)
ax2         = h2.subplots()
MoveFigure(h2,1250,10)

# Plot Initial Conditions
pl = Array{Any}(undef,nflds)
if initplot
  for j in 1:nflds
    if (plotfldi[j])
      pl[j] = ax2.plot(Geom.xm1[:],Q*fld[:,j],color=cm(j-1));
#      pl2   = ax2.plot(Geom.xm1[:],Q*fld[:,2],color=rgba1);
    end
  end  
  if ifphplot
    scat = ax1.plot(fld[:,1],fld[:,2],color="black") 
  end

  println("Press any key to continue")
  xin = readline()

# # Remove Plots  
#   for j in 1:nflds
#     if (plotfldi[j])
#       pl[j][1].remove();
#     end
#   end  
# #  pl[1].remove()
# #  pl2[1].remove()
#   if (ifphplot)
#     scat[1].remove()
#   end  
end

pause(0.01)






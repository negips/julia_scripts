
println("Newton Initialization")


# Fine grid interpolation operator
x2    = LinRange(Geom.xm1[1],Geom.xm1[end],500)
l2    = length(x2)
intpd = zeros(Float64,l2,lx1)
for i in 1:l2
  for j in 1:lx1
    intpd[i,j] = SpectralBases.Explaguerre(x2[i],j)
  end
end

agauss      = 0.0*Geom.xm1[:]
if ngauss == 1
  agauss    = exp.(-((Geom.xm1[:] .- x0)/σg).^2)
else  
  for i in 1:ngauss
    global agauss
    agauss  = agauss .+ ampgauss[i]*exp.(-((Geom.xm1[:] .- x0gauss[i])/σg).^2)
  end
end  

#agauss      = exp.(-((Geom.xm1[:] .- x0)/σg).^2)# .*(sign.(Geom.xm1[:] .- x0))
k0          = 7
asin        = sin.(2.0*π/(xe-xs)*k0*Geom.xm1[:])
acos        = cos.(2.0*π/(xe-xs)*k0*Geom.xm1[:])

ainit       = vimultg.*(QT*asin)
binit       = vimultg.*(QT*acos)

fld         = zeros(VT,ndof,nflds)

for j in 1:nflds
  fld[:,j]  = Amp0[j]*ainit .+ Off0[j] .+ σ0i[j]*(rand(ndof) .- 0.5)
end  


cm          = get_cmap("tab10");
rgba0       = cm(0); 
rgba1       = cm(1); 
rgba2       = cm(2);

if !@isdefined(h2)
  h2          = figure(num=2)
  ax2         = h2.subplots()
  MoveFigure(h2,1250,10)
else
  figure(h2)
  h2.clf()
  h2.show()
  ax2         = h2.subplots()
end

# Plot Initial Conditions
pl = Array{Any}(undef,nflds)
if initplot
  for j in 1:nflds
    if (plotfldi[j])
      pl[j] = ax2.plot(x2,intpd*fld[:,j],color=cm(j-1));
      h2.show()
#      pl2   = ax2.plot(Geom.xm1[:],Q*fld[:,2],color=rgba1);
    end
  end  
  if ifphplot
    scat = ax1.plot(fld[:,1],fld[:,2],color="black") 
    h1.show()
  end

  println("Press any key to continue")
  xin = readline()

end

pause(0.01)






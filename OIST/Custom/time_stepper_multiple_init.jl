
# Include the function files
include("sem_main.jl")
#include("Meinhardt.jl")
include("Dealias.jl")
include("../GetEXT.jl")
include("../GetBDF.jl")


println("Time-Stepping Initialization")

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
asin        = sin.(2.0*π/(xe-xs)*k0*Geom.xm1[:])
acos        = cos.(2.0*π/(xe-xs)*k0*Geom.xm1[:])

ainit       = vimultg.*(QT*agauss)
binit       = vimultg.*(QT*agauss)

nflds       = 2                                 # No of fields
fld         = zeros(VT,ndof,nflds)
dotfld      = zeros(VT,ndof,nflds)
fldlag      = zeros(VT,ndof,2,nflds)

Rhs         = zeros(VT,ndof,nflds)
Rhslag      = zeros(VT,ndof,2,nflds)

fld[:,1]    = ampB0*ainit .+ B0Off .+ σbi*(rand(ndof) .- 0.5)
fld[:,2]    = ampA0*binit .+ A0Off .+ σai*(rand(ndof) .- 0.5)

fldhist     = zeros(VT,npts,nsurf_save,nflds)
Thist       = zeros(VT,nsurf_save)

bdf         = zeros(Float64,4)
ext         = zeros(Float64,3)

cm          = get_cmap("tab10");
rgba0       = cm(0); 
rgba1       = cm(1); 
rgba2       = cm(2); 

time        = range(0.,step=dt,length=nsteps);

h2          = figure(num=2)
ax2         = h2.subplots()

t           = 0.

Thist[1] = t
for i in 1:nflds
  fldhist[:,1,i] = Q*fld[:,i]
end  

# Plot Initial Conditions
if initplot
  pl    = ax2.plot(Geom.xm1[:],Q*fld[:,1],color=rgba0);
  pl2   = ax2.plot(Geom.xm1[:],Q*fld[:,2],color=rgba1);
  if ifphplot
    scat = ax1.plot(fld[:,1],fld[:,2],color="black") 
  end
  println("Press any key to continue")
  xin = readline()
  pl[1].remove()
  pl2[1].remove()
  if (ifphplot)
    scat[1].remove()
  end  
end










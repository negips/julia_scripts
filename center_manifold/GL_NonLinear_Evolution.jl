println("Non-Linear evolution for Ginzburg Landau equations")

#include("Module_SEM1D/SEM1D.jl")
#using .SEM1D

#include("Module_StepperArnoldi/StepperArnoldi.jl")
##using .StepperArnoldi

#include("Module_CenterManifold/CenterManifold.jl")
##using .CenterManifold

using Peaks
using Statistics
using Random

include("$JULIACOMMON/RK4.jl")
include("NLGinzburgLandau.jl")
include("OP_RK4.jl")


lafs        = 16
lgfs        = 12

#include("GL_Setup.jl")
#-------------------------------------------------- 

close("all")

ifplot      = true 
verbose     = true
nsteps      = 500000
ifsave      = false
plotstep    = 20000
verbosestep = 20000
histstep    = 1000
nhist       = Int(nsteps/histstep)

rng   = MersenneTwister(1235)
vt    = Complex{Inp.Dtype}
zro   = vt(0)

#v     = randn(vt,ndof)*1.0;
v     = zeros(vt,ndof)
v    .= v .+ 0.01*exp.(-(xg.-5.75).^2)
v    .= v .+ conj.(v)

Hist  = zeros(vt,nhist)
Time  = zeros(Float64,nhist)

if (ifplot)
  hv  = figure(num=2,figsize=[8.,6.]);
  ax2 = gca()
  ax2.set_xlabel(L"x",fontsize=lafs)
  ax2.set_ylabel(L"A",fontsize=lafs)
end

cm    = get_cmap("tab10");
rgba0 = cm(0) 
rgba1 = cm(1) 
rgba2 = cm(2) 

dt    = 0.0001

NGL(x)= NLGinzburgLandau(OPg,ones(vt,ndof),x,δ[5],vt(0),vt(0),Inp.lbc,Inp.rbc)  

t     = Inp.Dtype(0)    # Time
i     = 0               # Istep

ifconv = false
println("Press x to stop. Any other key to continue...")
xin = readline()
if xin != "x"
  ifconv = false
else
  ifconv = true
end  

# Start iterations
println("Starting Iterations")

while (~ifconv)
  global v
  global t, i
  global plr,pli
  global OPg
  global hλ, hv, ax2

  local β

  β = one
  i = i + 1

  t = t + dt;

  # Apply BC
  SEM1D.SEM_SetBC!(v,Inp.lbc,Inp.rbc)
  # Non-linear Evolution
  v  = OP_RK4!(NGL,v,dt)

  # No Arnoldi iteration      
  if verbose && mod(i,verbosestep)==0
    println("Istep=$i, Time=$t")
  end

  if (mod(i,histstep) == 0)
    j = Int(i/histstep)
    Hist[j] = v[50]
    Time[j] = t
  end  

  # Plot the field  
  if (ifplot && mod(i,plotstep)==0)
    if (i>plotstep) 
      for lo in ax2.get_lines()
        lo.remove()
      end  
    end  
   
    pv1 = ax2.plot(xg,real.(v),linestyle="-",color=rgba0)
    pv2 = ax2.plot(xg,imag.(v),linestyle="--",color=rgba0)
    pv2 = ax2.plot(xg,abs.(v) ,linestyle="-",color=rgba1,linewidth=2)

    vmax = 1.2*maximum(abs.(v))
    vmin = -vmax
    dv   = abs(vmax-vmin)
    ax2.set_ylim((vmin,vmax))
    # ax2.set_ylim((-dv,dv))
    hv.show()    
  end 
  
  if i==nsteps
    break
  end  

end       # while ... 

if !ifconv
  h3  = figure(num=3,figsize=[8.,6.]);
  ax3 = gca()
  ax3.plot(Time,real.(Hist))
  ax3.set_xlabel(L"t",fontsize=lafs)
  ax3.set_ylabel(L"A",fontsize=lafs)
 
end


if (ifsave)
  # fname = "GL_BaseState5.jld2"
  # save(fname,"xg",xg,"basestate",v,"U",U,"γ",γ,"μx",μx,"μ0",μ0,"δ5",δ5,"δγ",δγ,"δμx",δμx);
  # println(fname*" saved.")
end 

println("Done.")












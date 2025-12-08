println("Non-Linear evolution for the Stuart Landau equations")

using Peaks
using Statistics
using Random

include("$JULIACOMMON/RK4.jl")
include("NLGinzburgLandau.jl")
include("OP_RK4.jl")


lafs        = 16
lgfs        = 12

# include("GL_Setup.jl")
#-------------------------------------------------- 

close("all")

cm          = get_cmap("tab10");
rgba0       = cm(0) 
rgba1       = cm(1) 
rgba2       = cm(2) 

ifplot      = false
histplot    = true
verbose     = true
nsteps      = 3000000
ifsave      = false
plotstep    = 20000
verbosestep = 10000
histstep    = 100
nhist       = Int(nsteps/histstep)
xhist       = true

dt          = 0.001


vt          = Complex{Inp.Dtype}
zro         = vt(0)

Hist        = zeros(vt,nhist,m)
Time        = zeros(Float64,nhist)
Mode_Ind    = [1]                         # Which mode to plot 
if xhist
  hist_x    = ForcingLocation()           # Location of history point
  hist_i    = argmin(abs.(xg .- hist_x))  # Index of history point

  Histx     = zeros(vt,nhist)
end  


# Stuart Landau
G1          = Khat
G2          = G_O2
G3          = G_O3
SL(x)       = StuartLandau3(G1,G2,G3,x)  

z           = zeros(vt,m)
rng         = Xoshiro(1235)
# Mode initial values
z[1]        = 1.0e-1*rand(rng,vt)
z[2]        = z[1]'
# Parameter Perturbations
z[n+1:n+p]  = zeros(vt,p)
# Harmonic Forcing Amplitude
z[n+p+1]    = 0.01*(1.0 + 0.0im)
z[n+p+2]    = z[n+p+1]'

ifconv = false
println("Press x to stop. Any other key to continue...")
xin = readline()
if xin == "x"
  nsteps = 0
end  


# Work Arrays
zwork = zeros(vt,m,5)

# Start iterations
t           = Inp.Dtype(0)    # Time
println("Starting Iterations")

for i in 1:nsteps
  global z
  global zwork
  global t
  global plr,pli
  global OPg
  global hλ, hv, ax2

  local β

  β = one

  t = t + dt;

  # Stuart Landau Evolution
  OP_RK4!(SL,z,dt,zwork)

  # Set conjugation correctly
  # z[2] = z[1]'

  # Print something  
  if verbose && mod(i,verbosestep)==0
    println("Istep=$i, Time=$t")
  end

  if (mod(i,histstep) == 0)
    j = Int(i/histstep)
    Hist[j,:]     = copy(z)
    Time[j]       = t

    # Get actual field value
    if (xhist)
      Histx[j] = Get_AsymptoticFieldx(hist_i,z,Vext,Y_O2,Y_O3)
    end  
  end  

  # # Plot the field  
  # if (ifplot && mod(i,plotstep)==0)
  #   if (i>plotstep) 
  #     for lo in ax2.get_lines()
  #       lo.remove()
  #     end  
  #   end  
  #  
  #   pv1 = ax2.plot(xg,real.(v),linestyle="-",color=rgba0)
  #   pv2 = ax2.plot(xg,imag.(v),linestyle="--",color=rgba0)
  #   pv2 = ax2.plot(xg,abs.(v) ,linestyle="-",color=rgba1,linewidth=2)
  # end 

end       # i in 1:nsteps 

if histplot && nsteps>0
  h3  = figure(num=3,figsize=[8.,6.]);
  ax3 = gca()
  ax3.plot(Time,real.(Hist[:,Mode_Ind]))
  ax3.set_xlabel(L"t",fontsize=lafs)
  ax3.set_ylabel(L"A",fontsize=lafs)

  if (xhist)
    h4  = figure(num=4,figsize=[8.,6.]);
    ax4 = gca()
    ax4.plot(Time,real.(Histx))
  end  
end


if (ifsave)
  # fname = "GL_BaseState5.jld2"
  # save(fname,"xg",xg,"basestate",v,"U",U,"γ",γ,"μx",μx,"μ0",μ0,"δ5",δ5,"δγ",δγ,"δμx",δμx);
  # println(fname*" saved.")
end 

println("Done.")












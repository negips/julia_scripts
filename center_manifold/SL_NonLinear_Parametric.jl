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
histstep    = 50
nhist       = Int(nsteps/histstep)
xhist       = true
vt          = Complex{Inp.Dtype}
zro         = vt(0)

dt          = 0.001
#θA          = [0.01; 0.1; 1.0]
θA          = [0.1]
nθ          = length(θA)

Hist_Mode   = zeros(vt,nhist,m,nθ)
Time        = zeros(Float64,nhist)
Mode_Ind    = [1]                         # Which mode to plot 

figsz        = [17.0, 6.0]

h3          = figure(num=3,figsize=figsz);
ax3         = gca()

if xhist
  hist_x    = 10.0                        # Location of history point
  hist_i    = argmin(abs.(xg .- hist_x))  # Index of history point

  Histx     = zeros(vt,nhist,nθ)

  h4        = figure(num=4,figsize=figsz);
  ax4       = gca()
end  


# Stuart Landau
G1          = Khat
G2          = G_O2
G3          = G_O3
SL(x)       = StuartLandau3(G1,G2,G3,x)  

println("Press x to stop. Any other key to continue...")
xin = readline()
if xin == "x"
  nsteps = 0
end  


for ik in 1:nθ

  # Initialize z
  θAmp        = θA[ik]

  z           = zeros(vt,m)
  rng         = Xoshiro(1235)
  # Mode initial values
  z[1]        = 1.0e-2*vt(1) #*rand(rng,vt)
  z[2]        = z[1]'
  # Parameter Perturbations
  z[n+1:n+p]  = zeros(vt,p)
  # Harmonic Forcing Amplitude
  z[n+p+1]    = θAmp*(1.0 + 0.0im)
  z[n+p+2]    = z[n+p+1]'
  
  # Work Arrays
  zwork       = zeros(vt,m,5)
  
  # Start iterations
  t           = Inp.Dtype(0)    # Time
  println("Starting Iterations")
  
  for i in 1:nsteps
  
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
      Hist_Mode[j,:,ik] = copy(z)
      Time[j]           = t
  
      # Get field value at point x = hist_x,
      # corresponding to array index hist_i
      if (xhist)
        Histx[j,ik] = Get_AsymptoticFieldx(hist_i,z,Vext,Y_O2,Y_O3)
      end  
    end  
  
  end       # i in 1:nsteps 

  if histplot && nsteps>0
    ax3.plot(Time,real.(Hist_Mode[:,Mode_Ind,ik]))
    ax3.set_xlabel(L"t",fontsize=lafs)
    ax3.set_ylabel(L"A",fontsize=lafs)
  
    if (xhist)
      leg = "|θ| = $(θAmp)"
      ax4.plot(Time,real.(Histx[:,ik]),label=leg)
    end  
  end

end         # ik in 1:nθ

ax4.set_xlim([2900.0,3000.0])
#ax4.set_ylim([-0.4,0.4])
ax4.legend(fontsize=lgfs,ncols=nθ)

println("Done.")












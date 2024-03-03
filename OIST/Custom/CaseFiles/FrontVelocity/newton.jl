println("Time-Stepping interface")

using PolynomialBases
using PyPlot,PyCall
using LinearAlgebra
using IterativeSolvers

include("$SRC/Dealias.jl")

include("newton_init.jl")

X = Geom.xm1[:];
X[end] = X[1]
QTX = QT*(X.*vimult)

# No dynamic phase here
ifdynplot         = false
ifplot            = iffldplot || ifphplot

C                 = 1.0       # Front Velocity
Rhs               = 0.0*fld

for i in 1:nsteps
  global fld,dotfld,Rhs,C
  global pl,pl2,scat,λpl
  global PlotContainers
  global framecount


  θpar      = θ0  # for G
  λpar      = λ0  # for F
  dotfld    = Flow(fld[:,1],fld[:,2],θpar,λpar)

  # Initial Residual
  for j in 1:1 # nflds
    rhs           = dotfld[:,j];
    lapfld        = Lg*fld[:,j];
    convfld       = Cg*fld[:,j];

    Rhs[:,j]      = lapfld .+ C*convfld .+ Bg.*rhs;

    # Boundary conditions: 0.0 at the right Boundary
    Rhs[end,j]    = 0.0
  end

  ϵ         = 0.001
  Cpert     = ϵ*rand()
  fldpert   = ϵ*rand(Float64,size(Fld[:,1]))
  fldnew    = fld[:,1] .+ fldpert

  dotfld2   = Flow(fldnew,fld[:,2],θpar,λpar)
  dotpert   = dotfld2[:,1] .- dotfld[:,1]

  # Action on perturbed system:
  rhspert   = Bg.*dotpert .+ Lg*fldpert .+ C*Cg*fldpert + Cpert*Cg*fld[:,1];







  if ifplot && mod(i,plotupd)==0
#   Remove old plots      
    if (i>plotupd)
       if (iffldplot)
         for j in 1:nflds
          if (plotfldi[j])
            pl[j][1].remove()
          end
        end  
       end  
       if (ifphplot)
         scat[1].remove()
       end
       if ifdynplot
         λpl[1].remove()
       end
    end

#   Add updated plots      
    if (iffldplot)
      for j in 1:nflds
        if (plotfldi[j])
          pl[j] = ax2.plot(Geom.xm1[:],Q*fld[:,j],color=cm(j-1));
        end
      end  
    end
#   Phase plot    
    if ifphplot
      scat = ax1.plot(fld[:,1],fld[:,2],color="black",linewidth=2) 
    end

    pause(0.001)
  end  

end










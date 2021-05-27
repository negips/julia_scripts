println("Main interface for 1D SEM")

using PolynomialBases
using PyPlot,PyCall
using LinearAlgebra
using IterativeSolvers


close("all")

# Include the function files
include("sem_main.jl")


# Temporal discretization
bdf1  = [ 1.  -1.  0.  0.]/1.;
bdf2  = [ 3.  -4.  1.  0.]/2.;
bdf3  = [11. -18.  9. -2.]/6.;
ex0   = [0. 1.  0. 0.];
ex1   = [0. 2. -1. 0.];
ex2   = [0. 3. -3. 1.];

# Local Matrices constructed in Sem_main.jl

Cg    = QT*Conv*Q    # Global Convection matrix
Lg    = QT*Lap*Q     # Global Laplacian matrix
Sg    = QT*Src*Q     # Global Src matrix
Fg    = QT*Fd*Q      # Global Feedback matrix
Bg    = QT*B         # Global Mass vector
Big   = 1.0./Bg      # Global inverse Mass vector

Oper  = similar(Bg)

xg    = QT*(vimult.*Geom.xm1[:])

Vinit = rand(Float64,ndof) + im*rand(Float64,ndof);
V     = copy(Vinit);
Vlag  = zeros(Complex,ndof,3);
Rlag  = zeros(Complex,ndof,2);

cm    = get_cmap("tab10");
rgba0 = cm(0); 
rgba1 = cm(1); 
rgba2 = cm(2); 

dt = 0.001;
plotupd = 100;

nsteps = 100000;

time = range(0.,step=dt,length=nsteps);

t = 0.
bc = zeros(Complex,1,ndof);
Rhs = similar(V)


for i in 1:nsteps
  global V,Vlag,Rlag
  global t
  global plr,pli
  global Oper,Rhs

  t = t + dt;

  println("Istep=$i, Time=$t")

  bdf = bdf3;
  ext = ex2;
     
  if i==1
    bdf = bdf1;
    ext = ex0;
  elseif i==2  
    bdf = bdf2;
    ext = ex1;
#  elseif i==3 
#    bdf = bdf3;
#    ext = ex2;
  end

# Build RHS
  rhs       = Cg*V            # Convection
  rhs       = rhs .+ Sg*V     # Source term
  rhs       = rhs .+ Fg*V     # Feedback

# Extrapolate Rhs  
  rhs1      = ext[2]*rhs + ext[3]*Rlag[:,1] + ext[4]*Rlag[:,2]

# Save Lag terms  
  Rlag[:,2] = Rlag[:,1]
  Rlag[:,1] = rhs

# Add contribution from lag terms  
  bdlag    = (1.0/dt).*Bg.*(bdf[2]*V + bdf[3]*Vlag[:,1] + bdf[4]*Vlag[:,2])
  Rhs      = rhs1 .- bdlag

  Vlag[:,3] = copy(Vlag[:,2])
  Vlag[:,2] = copy(Vlag[:,1])
  Vlag[:,1] = copy(V);

  Rhs[1]    = 0.       # Boundary condition at inlet

  Oper      = -Lg
  for j in 1:ndof
    Oper[j,j] = Oper[j,j] + bdf[1]/dt*Bg[j]
  end
  Oper[1,:] = bc
  Oper[1,1] = 1.0 + im*0.0

  soln         = gmres(Oper,Rhs);
  V            = copy(soln)

  if mod(i,plotupd)==0
    if (i>plotupd)
       plr[1].remove();
       pli[1].remove();
    end   
  
    plr = plot(xg,real.(V),color=rgba0);
    pli = plot(xg,imag.(V),color=rgba1,linestyle="--");

    pause(0.0001)
  end  

end
















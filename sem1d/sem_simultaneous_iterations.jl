println("Main interface for Time Stepper")

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

nkryl = 2

xg    = QT*(vimult.*Geom.xm1[:])

V     = rand(Float64,ndof,nkryl) + im*rand(Float64,ndof,nkryl);

# Orthogonalize

α           = sqrt(V[:,1]'*(Bg.*V[:,1]))
V[:,1]      = V[:,1]/α
for i in 2:nkryl
  global V
  local α
  h         = V[:,1:i-1]'*(Bg.*V[:,i])
  V[:,i]    = V[:,i] - V[:,1:i-1]*h
  α         = sqrt(V[:,i]'*(Bg.*V[:,i]))
  V[:,i]    = V[:,i]/α
end  

Vlag  = zeros(Complex,ndof,3,nkryl);
Rlag  = zeros(Complex,ndof,2,nkryl);

cm    = get_cmap("tab10");
rgba0 = cm(0) 
rgba1 = cm(1) 
rgba2 = cm(2) 

dt = 0.001
plotupd = 100
eigcal  = 500

nsteps = 100000

time = range(0.,step=dt,length=nsteps);

t = 0.
bc = zeros(Complex,1,ndof);
Rhs = similar(V[:,1])

verbose = false

for i in 1:nsteps
  global V,Vlag,Rlag
  global t
  global plr,pli

  t = t + dt;

  for ik in 1:nkryl
    global Oper,Rhs

    v = V[:,ik]  

    if verbose
      println("ik=$ik, Istep=$i, Time=$t")
    end  

    bdf = bdf3;
    ext = ex2;
       
    if i==1
      bdf = bdf1;
      ext = ex0;
    elseif i==2  
      bdf = bdf2;
      ext = ex1;
#    elseif i==3 
#      bdf = bdf3;
#      ext = ex2;
    end

#   Build RHS
    rhs       = Cg*v            # Convection
    rhs       = rhs .+ Sg*v     # Source term
    rhs       = rhs .+ Fg*v     # Feedback

#   Extrapolate Rhs  
    rhs1      = ext[2]*rhs + ext[3]*Rlag[:,1,ik] + ext[4]*Rlag[:,2,ik]

#   Save Lag terms  
    Rlag[:,2,ik] = Rlag[:,1,ik]
    Rlag[:,1,ik] = rhs

#   Add contribution from lag terms  
    bdlag    = (1.0/dt).*Bg.*(bdf[2]*v + bdf[3]*Vlag[:,1,ik] + bdf[4]*Vlag[:,2,ik])
    Rhs      = rhs1 .- bdlag

    Vlag[:,3,ik] = copy(Vlag[:,2,ik])
    Vlag[:,2,ik] = copy(Vlag[:,1,ik])
    Vlag[:,1,ik] = copy(v);

    Rhs[1]    = 0.       # Boundary condition at inlet

    Oper      = -Lg
    for j in 1:ndof
      Oper[j,j] = Oper[j,j] + bdf[1]/dt*Bg[j]
    end
    Oper[1,:] = bc
    Oper[1,1] = 1.0 + im*0.0

    soln         = gmres(Oper,Rhs);
    V[:,ik]      = copy(soln)

  end       # ik in 1:nkryl

  if mod(i,eigcal)==0
    println("Istep=$i, Time=$t")
   
    At = Vlag[:,1,:]'*V       # V'AV
    λ = eigvals(At)
    
    Lesshafft_λ = 1.0*im*λ
    display("$(Lesshafft_λ[1])")
    display("$(Lesshafft_λ[2])")
  end  

# Orthogonalization  
  β           = sqrt(V[:,1]'*(Bg.*V[:,1]))
  V[:,1]      = V[:,1]/β
  for i in 2:nkryl
    global V
    local β
    h         = V[:,1:i-1]'*(Bg.*V[:,i])
    V[:,i]    = V[:,i] - V[:,1:i-1]*h
    β         = sqrt(V[:,i]'*(Bg.*V[:,i]))
    V[:,i]    = V[:,i]/β
  end  

# Plot
  if mod(i,plotupd)==0
    if (i>plotupd)
       plr[1].remove();
       pli[1].remove();
    end   
  
    plr = plot(xg,real.(V[:,1]),color=rgba0);
    pli = plot(xg,imag.(V[:,1]),color=rgba1,linestyle="--");

    pause(0.0001)
  end
 
end















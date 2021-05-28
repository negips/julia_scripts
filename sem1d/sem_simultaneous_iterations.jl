println("Main interface for Time Stepper")

using PolynomialBases
using PyPlot,PyCall
using LinearAlgebra
using IterativeSolvers
using SpecialFunctions
using Roots

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

rcParams = PyPlot.PyDict(PyPlot.matplotlib."rcParams")

# Analytical Eigenvalues
ω1 = find_zero(airyai,(-3.0,-0.0))
ω2 = find_zero(airyai,(-5.0,-3.0))
ω3 = find_zero(airyai,(-6.0,-5.0))
ω4 = find_zero(airyai,(-7.0,-6.0))
ω5 = find_zero(airyai,(-8.0,-7.0))
ω6 = find_zero(airyai,(-9.5,-8.0))
ω7 = find_zero(airyai,(-10.5,-9.5))

ω  = [ω1, ω2, ω3, ω4, ω5, ω6, ω7]
U  = 6.0
γ  = 1.0 - im*1.0

Ω  = im*(U*U/8 .- U*U/4/γ .+ γ^(1.0/3.0)*(U^(4.0/3.0))/(160.0^(2.0/3.0))*ω)

rcParams["markers.fillstyle"] = "none"
hλ = figure(num=1,figsize=[8.,6.]);
ax1 = gca()
pΛ = plot(real.(Ω),imag.(Ω),linestyle="none",marker="o",markersize=8)


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

dt = 0.01
plotupd = Inf
eigcal  = 500

λn = zeros(Complex,nkryl)

nsteps = 10000000

time = range(0.,step=dt,length=nsteps);

t = 0.
bc = zeros(Complex,1,ndof);
Rhs = similar(V[:,1])

verbose = true
verbosestep = 100

rcParams["markers.fillstyle"] = "full"

for i in 1:nsteps
  global V,Vlag,Rlag
  global t
  global plr,pli

  t = t + dt;

  for ik in 1:nkryl
    global Oper,Rhs

    v = V[:,ik]  

    if verbose && mod(i,verbosestep)==0 && ik==1
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

# Calculate Eigenvalues of Reduced operator  
  if mod(i,eigcal)==0
    global hλ, ax1, λ, pλ
    if i==eigcal
#      hλ = figure(num=1,figsize=[8.,6.]);
#      ax1 = gca()
    else
#      ax1.clear()
      pλ[1].remove()
    end  

    println("Istep=$i, Time=$t")
   
    Ar = V'*(Cg .+ Sg .+ Fg .+ Lg)*V       # V'AV
    λ = eigvals(Ar)
    
    Lesshafft_λ = 1.0*im*λ
    for j in length(λ):-1:1  
      display("$(Lesshafft_λ[j])")
    end  

    pλ = ax1.plot(real.(Lesshafft_λ),imag.(Lesshafft_λ), linestyle="none",marker=".", markersize=8)
    pause(0.001)
  end  


# Plot
  if mod(i,plotupd)==0
    global hev, ax2
    if i==plotupd
      hev = figure(num=2,figsize=[8.,6.]);
      ax2 = gca()
    end  

    if (i>plotupd)
#       plr[1].remove();
#       pli[1].remove();
#      lo = ax2.get_lines()
      for lo in ax2.get_lines()
        lo.remove()
      end  
    end   

    for ik in 1:nkryl
      plr = ax2.plot(xg,real.(V[:,ik]),color=cm(ik-1));
      pli = ax2.plot(xg,imag.(V[:,ik]),color=cm(ik-1),linestyle="--");
    end  

    pause(0.0001)
  end
 
end















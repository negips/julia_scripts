println("Main interface for Time Stepper")

using PolynomialBases
using PyPlot,PyCall
using LinearAlgebra
using IterativeSolvers
using SpecialFunctions
using Roots

close("all")

# Include the function files
# include("sem_main.jl")


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

nkryl = 1

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

dt = 0.0001
plotupd = 5
eigcal  = 500

λn = zeros(Complex,nkryl)

nsteps = 10000000

time = range(0.,step=dt,length=nsteps);

t = 0.
bc = zeros(Complex,1,ndof);
Rhs = similar(V[:,1])

verbose = true
verbosestep = 100

reortho = 50

rcParams["markers.fillstyle"] = "full"

for i in 1:nsteps
  global V,Vlag,Rlag
  global t
  global plr,pli

  if i==1
    println("Starting time loop")
  end  

  t = t + dt;

  for ik in 1:nkryl
    global Opg

    v = V[:,ik]  

    if verbose && mod(i,verbosestep)==0 && ik==1
      println("ik=$ik, Istep=$i, Time=$t")
    end  

    OPg       = Cg .+ Sg .+ Lg
    for j in 1:ndof
      OPg[j,:] = OPg[j,:]./Bg[j]
    end  
    OPg[1,:] = bc
    OPg[1,1] = 1.0 + im*0.0        # Change operator for BC

#   Apply BC       
    v[1]      = 0.0 + im*0.0
   
#   k1
    v1        = v .+ dt/2.0*OPg*v
    v2        = v .+ dt/2.0*OPg*v1
    v3        = v .+ dt*OPg*v2
    v4        = v .+ dt/6.0*(OPg*(v .+ 2.0*v1 .+ 2.0*v2 .+ v3))

    v         = v4

    V[:,ik]      = copy(v)

  end       # ik in 1:nkryl

# Orthogonalization && eig calculation
  if mod(i,reortho)==0
#   Calculate Eigenvalues of Reduced operator  
    global hλ, ax1, λ, pλ, Ar
   
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

    if i==reortho
#      hλ = figure(num=1,figsize=[8.,6.]);
#      ax1 = gca()
    else
#      ax1.clear()
      pλ[1].remove()
    end  
   
    Ar = V'*(Cg .+ Sg .+ Fg .+ Lg)*V       # V'AV
    λ = eigvals(Ar)
    
    Lesshafft_λ = 1.0*im*λ
    for j in length(λ):-1:1  
      display("$(Lesshafft_λ[j])")
    end  

    pλ = ax1.plot(real.(Lesshafft_λ),imag.(Lesshafft_λ), linestyle="none",marker=".", markersize=8)

    fac = 0.2
    o1 = minimum(real.(Ω))
    l1 = minimum(real.(Lesshafft_λ))
    x1 = min(l1,o1)

    o2 = maximum(real.(Ω))
    l2 = maximum(real.(Lesshafft_λ))
    x2 = max(l2,o2)

    dx = (x2-x1)*fac
    x1 = x1 - dx
    x2 = x2 + dx

    o1 = minimum(imag.(Ω))
    l1 = minimum(imag.(Lesshafft_λ))
    y1 = min(l1,o1) 

    o2 = maximum(imag.(Ω))
    l2 = maximum(imag.(Lesshafft_λ))
    y2 = max(l2,o2)

    dy = (y2-y1)*fac
    y1 = y1 - dy
    y2 = y2 + dy


    ax1.set_xlim((x1,x2))  
    ax1.set_ylim((y1,y2))  

   
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















# Include the function files
function StepArn(L::AbstractMatrix{T},B::AbstractVector{S},StpInp::StepperInput,ArnInp::ArnoldiInput,lbc::Bool,rbc::Bool) where {S,T<:Number}

  dtype           = T

  V,H             = ArnKrylovInit(StpInp,ArnInp;Dtype=T)
  v0              = ArnInitVector(ArnInp.vlen,lbc,rbc,Dtype=T) 

  nev             = ArnInp.nev
  eigshift        = ArnInp.eigshift
  tol             = ArnInp.tol
  ngs             = ArnInp.ngs

  ifarnoldi       = ArnInp.ifarnoldi
  verbose         = ArnInp.ifverbose
  ifeigshift      = ArnInp.ifeigshift
  ifadjoint       = StpInp.ifadjoint         
  nsteps          = StpInp.nsteps            # Stepper Phase steps
  dt              = StpInp.timestep          # Time step

  nkryl           = 0
  lkryl           = ArnInp.lkryl
  block           = ArnInp.bsize 
  h,θ,v           = ArnUpd(V,block,B,v0,nkryl,ngs)
  V[:,1]          = v
  nkryl           = 0

  Rhs             = similar(v)

  ifconv          = false
  t               = T(0)*dt       # Time
  i               = 0             # Istep

  maxouter_it     = ArnInp.outer_iterations
  major_it        = 1

  v1              = zeros(T,ArnInp.vlen) 
  v2              = zeros(T,ArnInp.vlen) 
  v3              = zeros(T,ArnInp.vlen)
  Bi              = 1.0./B          # Inverse Mass (Vector)

  # Eigenvalue Shift
  if ifeigshift
    Ω = exp(eigshift*nsteps*dt)
    println("EigShift: $(eigshift), Ω= $Ω, |Ω| = $(abs(Ω))")
  else
    Ω = 0.0
  end

  # Start iterations
  ifdirect = true
  println("Starting Stepper/Arnoldi Iterations")
  while (~ifconv)

    for i in 1:nsteps
      t = t + dt;
      # Apply BC
      StpArn_SetBC!(v,lbc,rbc)
      # RK4!(v,L,v1,v2,v3,dt)
      BiRK4!(v,L,Bi,v1,v2,v3,dt)
    end  
  
    # Expand Krylov space
    V,H,nkryl,β,major_it = IRAM!(V,H,B,v,nkryl,lkryl,major_it,nev,Ω,ngs)
    v   = V[:,nkryl]
  
    if (major_it>maxouter_it)
      break
    end  
  
    if (verbose)
      @printf "Major Iteration: %3i/%3i, Krylov Size: %3i/%3i, β: %12e\n" major_it maxouter_it nkryl lkryl β
    end
    if (β < tol)
      @printf "Stopping Iteration, β: %12e\n" β
      ifconv = true
    end  
  end      # while ... 


  # Plotting after convergence
  Hr  = H[1:nev,1:nev]
  F   = eigen(Hr)
  evs = F.values
  DT  = dt*nsteps
 

  ritz = evs
  
  λr   = log.(abs.(evs))/DT
  λi   = atan.(imag(evs),real.(evs))/DT
  
  λ  = λr .+ im*λi
  # Eigenvectors  
  eigvec = V[:,1:nev]*F.vectors

  # Structured output
  arnout = ArnoldiOutput(λ,eigvec,ritz,ifconv,tol) 
  
  return arnout
end
#---------------------------------------------------------------------- 
function RestrictedStepArn(L::AbstractMatrix{T},B::AbstractVector{S},Vr::AbstractMatrix{T},Wr::AbstractMatrix{T},StpInp::StepperInput,ArnInp::ArnoldiInput,lbc::Bool,rbc::Bool) where {S,T}

  dtype           = T

  V,H             = ArnKrylovInit(StpInp,ArnInp;Dtype=T)
  v0              = ArnInitVector(ArnInp.vlen,lbc,rbc,Dtype=T)

  # Remove components
  ObliqueSubspaceRemoval!(v0,Vr,Wr,B,ArnInp.ngs) 
  StpArn_SetBC!(v0,lbc,rbc)

  nev             = ArnInp.nev
  tol             = ArnInp.tol
  eigshift        = ArnInp.eigshift
  ngs             = ArnInp.ngs
  tol             = ArnInp.tol

  ifarnoldi       = ArnInp.ifarnoldi
  verbose         = ArnInp.ifverbose
  ifeigshift      = ArnInp.ifeigshift
  ifadjoint       = StpInp.ifadjoint         
  nsteps          = StpInp.nsteps            # Stepper Phase steps
  dt              = StpInp.timestep          # Time step

  nkryl           = 0
  lkryl           = ArnInp.lkryl
  block           = ArnInp.bsize 
  h,θ,v           = ArnUpd(V,block,B,v0,nkryl,ngs)
  V[:,1]          = v
  nkryl           = 0

  Rhs             = similar(v)

  ifconv          = false
  t               = T(0)*dt       # Time
  i               = 0             # Istep

  maxouter_it     = ArnInp.outer_iterations
  major_it        = 1

  v1              = zeros(T,ArnInp.vlen) 
  v2              = zeros(T,ArnInp.vlen) 
  v3              = zeros(T,ArnInp.vlen)
  Bi              = 1.0./B        # Inverse Mass (Vector)

  # Eigenvalue Shift
  if ifeigshift
    Ω = exp(eigshift*nsteps*dt)
    println("EigShift: $(eigshift), Ω= $Ω, |Ω| = $(abs(Ω))")
  else
    Ω = 0.0
  end

  # Start iterations
  ifdirect = true
  println("Starting Stepper/Arnoldi Iterations")
  while (~ifconv)

    for i in 1:nsteps
      t = t + dt;
      # Apply BC
      # StpArn_SetBC!(v,lbc,rbc)
      RestrictedBRK4!(v,L,B,Vr,Wr,lbc,rbc,v1,v2,v3,dt)
    end  
  
    # Expand Krylov space
    V,H,nkryl,β,major_it = IRAM!(V,H,B,v,nkryl,lkryl,major_it,nev,Ω,ngs)
    v   = V[:,nkryl]
  
    if (major_it>maxouter_it)
      break
    end  
  
    if (verbose)
      @printf "Major Iteration: %3i/%3i, Krylov Size: %3i/%3i, β: %12e\n" major_it maxouter_it nkryl lkryl β
    end
    if (β < tol)
      @printf "Stopping Iteration, β: %12e\n" β
      ifconv = true
    end  
  end      # while ... 


  # Plotting after convergence
  Hr  = H[1:nev,1:nev]
  F   = eigen(Hr)
  evs = F.values
  DT  = dt*nsteps

  ritz = evs
  λr   = log.(abs.(evs))/DT
  λi   = atan.(imag(evs),real.(evs))/DT
  
  λ  = λr .+ im*λi
  # Eigenvectors  
  eigvec = V[:,1:nev]*F.vectors

  # Structured output
  arnout = ArnoldiOutput(λ,eigvec,ritz,ifconv,tol) 
  
  return arnout
end
#---------------------------------------------------------------------- 






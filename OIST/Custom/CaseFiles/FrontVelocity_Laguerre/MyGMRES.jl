#---------------------------------------------------------------------- 
function MyGMRES!(sol::AbstractArray{T},r::AbstractArray{T},F::AbstractMatrix{T},V::AbstractArray{T},tol::Float64,maxoit::Integer) where {T <: Number}

  l,c    = size(V)
  nr     = ndims(r)
  ns     = ndims(sol)
  n      = length(r)

  @assert nr == 1 "Implementation for only single rhs"
  @assert ns == 1 "Implementation for only single rhs"

  @assert n  == l "Unequal Krylov space and input residual lengths"

  lk  = c-1
  global Hes = zeros(T,lk,lk)
  global Giv = Vector{LinearAlgebra.Givens}(undef,lk)
  global rhs = zeros(T,lk+1)     # Residuals
  α   = zeros(T,lk+1)

  zro = T(0)

  wk  = zeros(T,l)        # Work array. Maybe I should pass it?
  r2  = zeros(T,l)        # Work array. Maybe I should pass it?

  v   = zeros(T,n)

  res = norm(r)
 
  # Should remove this if we want to start with an initial guess
  fill!(sol,zro)

  oit  = 0
  while oit<maxoit && res>tol

    # Outer Iterations
    oit     = oit + 1 

    wk      = F*sol
    r2      = r .- wk

    res     = norm(r2)

    if res == 0.0
      println("Exact Initial Guess. Exitting in MyGMRES!()")
      return nothing
    end

    # if oit == 1
      println("Initial Residual: $res")
    # end

    # Initialize rhs
    fill!(rhs,zro)
    rhs[1]  = res

    wk      = r2/res
    copyto!(V,1,wk,1,n)

    nk = 0
    for k in 1:lk
      nk     += 1 
      v       = V[:,k]
      wk      = F*v          # Using Matrix
      α[1:k]  = V[:,1:k]'*wk
      wk      = wk .- V[:,1:k]*α[1:k]
      β       = norm(wk)

      α[k+1]  = β

      # Apply old Givens Rotations to new vector
      for j in 1:k-1
        α     = Giv[j]*α
      end  

      # Calculate new Givens Rotation
      Giv[k],γ  = givens(α,k,k+1)
      α[k]      = γ
      j         = (k-1)*(lk)+1
      copyto!(Hes,j,α,1,k)

      # Apply Givens to the rhs
      rhs     = Giv[k]*rhs

      res     = abs(rhs[k+1])
      if res<tol
        break
      end

      # Update Krylov space
      wk      = wk/β
      j       = (k*l)+1
      copyto!(V,j,wk,1,n)
    end     # for k in 1:lk

    # Construct Solution
    fill!(α,zro)
    for k in nk:-1:1
      temp = rhs[k]
      for j in k+1:nk
        temp = temp - Hes[k,j]*α[j]
      end
      α[k] = temp/Hes[k,k]
    end 

    wk   = V[:,1:nk]*α[1:nk]

    for i in 1:n
      sol[i] = sol[i] + wk[i]
    end  

  end       # while oit<maxoit && res>tol

  wk      = F*sol
  r2      = r .- wk
  res     = norm(r2)
  println("Final Residual: $res")

  return nothing
end
#---------------------------------------------------------------------- 
function MyGMRES!(sol::AbstractArray{T},r::AbstractArray{T},F::Function,V::AbstractArray{T},tol::Float64,maxoit::Integer) where {T <: Number}

  l,c    = size(V)
  nr     = ndims(r)
  ns     = ndims(sol)
  n      = length(r)

  @assert nr == 1 "Implementation for only single rhs"
  @assert ns == 1 "Implementation for only single rhs"

  @assert n  == l "Unequal Krylov space and input residual lengths"

  lk  = c-1
  Hes = zeros(T,lk,lk)
  Giv = Vector{LinearAlgebra.Givens}(undef,lk)
  rhs = zeros(T,lk+1)     # Residuals
  α   = zeros(T,lk+1)

  zro = T(0)

  wk  = zeros(T,l)        # Work array. Maybe I should pass it?
  r2  = zeros(T,l)        # Work array. Maybe I should pass it?

  v   = zeros(T,n)

  res = norm(r)
 
  # Should remove this if we want to start with an initial guess
  fill!(sol,zro)

  oit  = 0
  while oit<maxoit && res>tol

    # Outer Iterations
    oit     = oit + 1 

    wk      = F(sol)
    r2      = r .- wk

    res     = norm(r2)

    if res == 0.0
      println("Exact Initial Guess. Exitting in MyGMRES!()")
      return nothing
    end

    if oit == 1
      println("Initial Residual: $res")
    end

    # Initialize rhs
    fill!(rhs,zro)
    rhs[1]  = res

    wk      = r2/res
    copyto!(V,1,wk,1,n)

    nk = 0
    for k in 1:lk
      nk     += 1 
      v       = V[:,k]
      wk      = F(v)          # Using Anonymous function
      α[1:k]  = V[:,1:k]'*wk
      wk      = wk .- V[:,1:k]*α[1:k]
      β       = norm(wk)

      α[k+1]  = β

      # Apply old Givens Rotations to new vector
      for j in 1:k-1
        α     = Giv[j]*α
      end  

      # Calculate new Givens Rotation
      Giv[k],γ  = givens(α,k,k+1)
      α[k]      = γ
      j         = (k-1)*(lk)+1
      copyto!(Hes,j,α,1,k)

      # Apply Givens to the rhs
      rhs     = Giv[k]*rhs

      res     = abs(rhs[k+1])
      if res<tol
        break
      end

      # Update Krylov space
      wk      = wk/β
      j       = (k*l)+1
      copyto!(V,j,wk,1,n)
    end     # for k in 1:lk

    # Construct Solution
    fill!(α,zro)
    for k in nk:-1:1
      temp = rhs[k]
      for j in k+1:nk
        temp = temp - Hes[k,j]*α[j]
      end
      α[k] = temp/Hes[k,k]
    end 

    wk   = V[:,1:nk]*α[1:nk]

    for i in 1:n
      sol[i] = sol[i] + wk[i]
    end  

  end       # while oit<maxoit && res>tol

  wk      = F(sol)
  r2      = r .- wk
  res     = norm(r2)
  println("Final Residual: $res")

  return nothing
end
#---------------------------------------------------------------------- 











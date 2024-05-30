function SetZero!(A::AbstractMatrix{T},tol::Float64) where {T}

  zro = T(0)
  for i in eachindex(A)
    if abs(A[i])<tol
      A[i] = zro
    end
  end

  return nothing
end  
#---------------------------------------------------------------------- 

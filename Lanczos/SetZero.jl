function SetZero!(A,tol)

  n = length(A)
  el = eltype(A[1])
  zro = el(0)
  for i in 1:n
    if abs(A[i])<tol
      A[i] = zro
    end
  end

end  
#---------------------------------------------------------------------- 

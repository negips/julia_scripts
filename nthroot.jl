function nthroot(A::Float64,n::Int,tol::Float64)

  xk = 2.0
  xn = xk^n
  k  = 0
  kmax = 100
  relerr = (xn/A-1.0)

  while abs(relerr)>tol

    k += 1
    if k>kmax
      println("Unconverged nth root!", kmax)
      break
    end

    xk = xk*(n-1.0)/n + (A/n)*(xk/xn)
    xn = xk^n
    relerr = (xn/A-1.0)
   
  end

  return xk
end  

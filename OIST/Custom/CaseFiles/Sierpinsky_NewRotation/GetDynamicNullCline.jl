function GetDynamicNullCline(f,y,λ)

  # Assuming f is linear in x

  f00       = f(0.0,0.0,λ)          # inhomogeneous
  cx        = f(1.0,0.0,λ) - f00    # x coefficient
  xnull     = -f(0.0,y,λ)./cx

  return xnull
end



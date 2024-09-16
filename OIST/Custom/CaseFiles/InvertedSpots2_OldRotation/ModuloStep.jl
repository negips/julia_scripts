function ModuloStep(x,L::Float64)

  xfloor = floor.(x./L)
  xmod   = x/L .- xfloor .- 0.5

  step = sign.(xmod)
  
  return step
end  

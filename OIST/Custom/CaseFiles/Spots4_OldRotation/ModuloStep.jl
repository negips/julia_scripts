function ModuloStep(x,L::Float64)

  L2 = 2.0*L

  xfloor = floor.(x./L2)
  xmod   = x/L2 .- xfloor .- 0.5

  step = sign.(xmod)
  
  return step
end  

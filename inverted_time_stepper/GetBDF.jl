#!/usr/bin/jl

function GetBDF!(bdfk,ord::Int)
 
  # Temporal discretization
  if ord==1
    bdfk[1] =  1.0 
    bdfk[2] = -1.0
    bdfk[3] =  0.0
    bdfk[4] =  0.0
   elseif ord==2
    bdfk[1] =  3.0/2.0
    bdfk[2] = -4.0/2.0 
    bdfk[3] =  1.0/2.0
    bdfk[4] =  0.0
  elseif ord==3
    bdfk[1] =  11.0/6.0
    bdfk[2] = -18.0/6.0 
    bdfk[3] =  9.0/6.0
    bdfk[4] = -2.0/6.0
  else
    @assert false "Order greater than 3 not implemented."
  end

  return
end  
    


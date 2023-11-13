#!/usr/bin/jl

function GetEXT!(extk,ord::Int)
 
  # Temporal discretization
  if ord==1
    extk[1] =  1.0 
    extk[2] =  0.0 
    extk[3] =  0.0
   elseif ord==2
    extk[1] =  2.0
    extk[2] = -1.0 
    extk[3] =  0.0 
  elseif ord==3
    extk[1] =  3.0 
    extk[2] = -3.0 
    extk[3] =  1.0
  else
    @assert false "Order greater than 3 not implemented."
  end

  return
end  
    


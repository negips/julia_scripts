#!/usr/bin/julia
@doc raw"""
      function GetBDF(Ord::Int)

      Backward Difference Coefficients for first order time derivative

"""
function GetBDF(ord::Int)
 
  bdfk = zeros(Float64,4)
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

  return bdfk
end  

#---------------------------------------------------------------------- 
@doc raw"""
      function GetBDF2(Ord::Int)

      Backward Difference Coefficients for second order time derivative

"""
function GetBDF2(ord::Int)
 
  @assert ord<3 "Order greater than 2 not implemented."
  @assert ord>0 "Invalid Order $ord."

  bdfk = zeros(Float64,4)
  # Temporal discretization
  if ord==1
    bdfk[1] =  1.0 
    bdfk[2] = -2.0
    bdfk[3] =  1.0
    bdfk[4] =  0.0
  elseif ord==2
    bdfk[1] =  2.0
    bdfk[2] = -5.0
    bdfk[3] =  4.0
    bdfk[4] = -1.0
  end

  return bdfk
end  
#---------------------------------------------------------------------- 
@doc raw"""
      function GetEXT(Ord::Int)

      Extrapolation Coefficients 

"""
function GetEXT(ord::Int)
 
  # Extrapolation Coefficients
  extk = zeros(Float64,3)
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

  return extk
end
#---------------------------------------------------------------------- 



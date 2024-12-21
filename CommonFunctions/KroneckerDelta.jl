# Kronecker Delta function
function KroneckerDelta(a::Int,b::Int)
  if a == b
    return 1
  else
    return 0
  end
end  

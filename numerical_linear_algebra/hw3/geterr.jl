using LinearAlgebra
function geterr(A)
# Input:
# A - input Matrix

  local lA
  local err
  n = size(A);
  if n[1] != n[2]
    println("A is not a square matrix $n")
  end
  
  lA = copy(LowerTriangular(A));
  for i in 1:n[1]
    lA[i,i]= 0;
  end
  err=maximum(maximum(abs.(lA)));

  return err

end


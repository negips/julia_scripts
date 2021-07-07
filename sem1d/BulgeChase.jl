# Bulge Chase Algorithm
function BulgeChase(H::Matrix,μ::Vector,nμ::Int)

  # nμ      - No of Shifts
  # μ       - Shifts
  # H       - Matrix

  x0    = 0. *H[:,1]
  x0[1] = 1.0
  type  = typeof(H[1,1])
  x     = zeros(type,length(x0)) 

  # calculate x = p(A)*e1 = (A - μnI)...(A-μ2I)(A-μ1I)*e1
  for i in 1:nμ
    global x0, x
    for j in 1:i+1
      x[j] = H[j,1:i]*x0[1:i] - μ[i]*x0[j]
    end  
    x0 = x
  end  

  return x0
end  

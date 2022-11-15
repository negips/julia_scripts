function UpperHessenbergReduction(A::Matrix)

# Only for Square matrices 
  r,c = size(A)
  el  = eltype(A[1])
  zro = el(0)

  V   = Matrix{el}(I,r,c)
  B   = zeros(el,r,c)

  if r!=c
    display("H needs to be a square matrix: size(A)= $r,$c")
  end

  v0 = zeros(el,r)
  for i in c-2:-1:c-2 #c-2
    x             = copy(A[:,i])
    v             = HouseHolderReflector(x,i+1) 
    Q             = (I - el(2)*v*v')
    mul!(B,Q,A)
#    A            .= B
    mul!(A,B,Q')
    V[:,i]        = v 
  end

  return V
end
#----------------------------------------------------------------------
function HouseHolderReflector(x,i)

  el              = eltype(x[1])
  r               = length(x)
  v               = zeros(el,r)
  v[i]            = el(1)
  v[i:r]         .= sign(x[i])*norm(x[i:r])*v[i:r] .+ x[i:r]
  v               = v./norm(v)

  return v
end
#---------------------------------------------------------------------- 

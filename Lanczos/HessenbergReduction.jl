function UpperHessenbergReduction(A::Matrix)

# Only for Square matrices 
  r,c = size(A)
  el  = eltype(A[1])
  zro = el(0)

  V   = Matrix{el}(I,r,c)
  B   = zeros(el,r,c)

  if r!=c
    display("A needs to be a square matrix: size(A)= $r,$c")
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
function LowerHessenbergtoTriDiagonal!(H::Matrix)

  el = eltype(H[1])
  one = el(1)

  tol = 1000*eps(abs.(one))

  r,c = size(H)
  if r!=c
    display("H needs to be a square matrix: size(H)= $r,$c")
  end

  W  = zeros(el,r,c)
  V  = zeros(el,r,c)

  for i in 1:c-2
    w       = transpose(H[i,:])
    w[1:i]  = zeros(el,1,i)
    v       = H[:,i]
    v[1:i]  = zeros(el,i)
    v[i+1]  = one/w[i+1]

    if abs(w[i+1])<tol
      wa    = abs(w[i+1])
      println("Vanishing row abs(w)=$wa")
    end  
   
    β       = (w*H[:,i])[1]
    if abs(β)<tol
      βa    = abs(β)
      println("Orthogonal vectors: abs(β)=$βa")
    end  
    w       = w/β
    Q       = I - v*w
    A       = Q*H
    V[:,i+1] = -v*w[i+1]
    V[i+1,i+1] = one + V[i+1,i+1]

#   Build and store inverse Q as vectors in W
    W[i+1,i+1]  = one/V[i+1,i+1]
    for j in i+2:r
      W[j,i+1] = - V[j,i+1]/V[i+1,i+1] 
    end

    Q[:,i+1]  = W[:,i+1]
    H        .= A*Q
    
  end

  return V,W
end
#----------------------------------------------------------------------
function UpperHessenbergtoTriDiagonal!(H::Matrix)

  el = eltype(H[1])
  one = el(1)

  tol = 1000*eps(abs.(one))

  r,c = size(H)
  if r!=c
    display("H needs to be a square matrix: size(H)= $r,$c")
  end

  W  = zeros(el,r,c)
  V  = zeros(el,r,c)

  for i in 1:r-2
    w       = transpose(H[i,:])
    w[1:i]  = zeros(el,1,i)
    v       = H[:,i]
    v[1:i]  = zeros(el,i)
    w[i+1]  = one/v[i+1]

    if abs(v[i+1])<tol
      va    = abs(v[i+1])
      println("Vanishing column abs(v)=$va")
    end  

    β       = (transpose(H[i,:])*v)[1]
    if abs(β)<tol
      βa    = abs(β)
      println("Orthogonal vectors: abs(β)=$βa")
    end  
   
    v       = v/β
    Q       = I - v*w
    A       = H*Q
    W[i+1,:] = -v[i+1]*w
    W[i+1,i+1] = one + W[i+1,i+1]

#    H      .= A 

#   Build and store inverse Q as vectors in W
    V[i+1,i+1]  = one/W[i+1,i+1]
    for j in i+2:c
      V[i+1,j] = - W[i+1,j]/W[i+1,i+1] 
    end

    Q[i+1,:]  = V[i+1,:]
    H      .= Q*A
    
  end

  return V,W
end
#----------------------------------------------------------------------






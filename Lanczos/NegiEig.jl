# Bulge Chase Algorithm with Oblique projectors
function NegiAlg(T::Matrix,λ)
# Also known as BulgeChase Algorithm  

  # H       - Tridiagonal Matrix
  # λ       - Shift

  r,c   = size(T)

# Create a bulge in the matrix
  T,v,w     = CreateLowerBulgeOblique(T,λ)
  Vi,Wi     = LowerHessenbergtoTriDiagonal!(T)

# Collect Left multipliers 
  V = Vi
# Collect Right multipliers
  W = Wi

  return T,V,W
end
#----------------------------------------------------------------------
# Lower Bulge Chase Algorithm with Oblique projectors
function NegiAlg2(T::Matrix,λ)
# Also known as BulgeChase Algorithm  

  # H       - Tridiagonal Matrix
  # λ       - Shift

  r,c   = size(T)

# Create a bulge in the matrix
  T,v,w     = CreateLowerBulgeOblique(T,λ)
  Vi,Wi     = ChaseBulgeTriDiagonal2!(T)

# Collect Left multipliers 
  V = Vi
# Collect Right multipliers
  W = Wi

  return T,V,W
end
#----------------------------------------------------------------------
# Upper Bulge Chase Algorithm with Oblique projectors
function NegiAlg3(T::Matrix,λ)
# Also known as BulgeChase Algorithm  

  # H       - Tridiagonal Matrix
  # λ       - Shift

  r,c   = size(T)

# Create a bulge in the matrix
  T,v,w     = CreateUpperBulgeOblique(T,λ)
  Vi,Wi     = ChaseBulgeTriDiagonal3!(T)

# Collect Left multipliers 
  V = Vi
# Collect Right multipliers
  W = Wi

  return T,V,W
end
#----------------------------------------------------------------------

function CreateLowerBulgeOblique(T::Matrix,λ)

  # T       - Tri-diagonal Matrix
  # λ       - Shifts

  T         = T - λ*I
  el        = eltype(T[1])
  zro       = el(0)
  one       = el(1)
  r,c       = size(T)
  x         = zeros(el,r)
  x[1]      = zro
  x[2]      = T[2,1]
  x         = x
  y         = transpose(T[1,:])
  α         = (y*x)[1]
  y         = y/α      # This makes <y,x> = 1.0
  β         = (y*T[:,1])[1]
  y         = y/β      # => <y,x> = 1.0/β
  V         = I - (x*y)
  Vi        = I + (x*y)/(one - one/β)

# This creates the Bulge  
# A = Q0*H*Q0  
  tmp       = V*T
  T        .= tmp*Vi
  T         = T + λ*I

  return T,x,y
end

#----------------------------------------------------------------------
function CreateUpperBulgeOblique(T::Matrix,λ)

  # T       - Tri-diagonal Matrix
  # λ       - Shifts

  T         = T - λ*I
  el        = eltype(T[1])
  zro       = el(0)
  one       = el(1)
  r,c       = size(T)
  y         = zeros(el,1,c)
  y[1]      = zro
  y[2]      = T[1,2]
  x         = T[:,1]
  α         = (y*x)[1]
  x         = x/α      # This makes <y,x> = 1.0
  β         = (transpose(T[1,:])*x)[1]
  x         = x/β      # => <y,x> = 1.0/β
  W         = I - (x*y)
  Wi        = I + (x*y)/(one - one/β)

# This creates the Bulge  
# A = Q0*H*Q0  
  tmp       = T*W
  T        .= Wi*tmp
  T         = T + λ*I

  return T,Wi,W
end

#----------------------------------------------------------------------
function LowerHessenbergtoTriDiagonal2!(H::Matrix)

  el = eltype(H[1])
  one = el(1)

  tol = 1000*eps(abs.(one))

  r,c = size(H)
  if r!=c
    display("H needs to be a square matrix: size(H)= $r,$c")
  end

  W  = Matrix{vt}(I,r,c)
  V  = Matrix{vt}(I,r,c)

  one = vt(1)

  for i in 1:c-2
    w       = transpose(H[i,:])
    w[1:i]  = zeros(el,1,i)
    v       = H[:,i]
    v[1:i]  = zeros(el,i)
    v[i+1]  = one/w[i+1]      # <w,v> = 1.0

    if abs(w[i+1])<tol
      wa    = abs(w[i+1])
      println("Vanishing row abs(w)=$wa")
    end  
   
    β       = (w*H[:,i])[1]
    if abs(β)<tol
      βa    = abs(β)
      println("Orthogonal vectors: abs(β)=$βa")
    end  
    w       = w/β             # <w,v> = 1.0/β

#   Left Multiplication      
    Q       = I - v*w
    A       = Q*H
    W       = Q*W

    β       = (w*v)[1]
#   Right multiplication    
    Q       = I + (v*w)/(one - β)
    H      .= A*Q
    V       = V*Q
    
  end

  return V,W
end
#----------------------------------------------------------------------
function UpperHessenbergtoTriDiagonal2!(H::Matrix)

  el = eltype(H[1])
  one = el(1)

  tol = 1000*eps(abs.(one))

  r,c = size(H)
  if r!=c
    display("H needs to be a square matrix: size(H)= $r,$c")
  end

  W  = Matrix{vt}(I,r,c)
  V  = Matrix{vt}(I,r,c)

  ngs = 2

  one = vt(1)

  for i in 1:c-2
    w       = transpose(H[i,:])
    w[1:i]  = zeros(el,1,i)
    v       = H[:,i]
    v[1:i]  = zeros(el,i)
    v[i+1]  = one/w[i+1]      # <w,v> = 1.0

    if abs(v[i+1])<tol
      va    = abs(v[i+1])
      println("Vanishing row column(v)=$va")
    end  
   
    β       = (transpose(H[i,:])*v)[1]
    if abs(β)<tol
      βa    = abs(β)
      println("Orthogonal vectors: abs(β)=$βa")
    end  

    v       = v/β             # <w,v> = 1.0/β
#   Right Multiplication      
    Q       = I - v*w
    A       = H*Q
    V       = V*Q

    β       = (w*v)[1]
#   Left multiplication    
    Q       = I + (v*w)/(one - β)
    H      .= Q*A
    W       = Q*W
  end

  return V,W
end
#----------------------------------------------------------------------

function ChaseBulgeTriDiagonal!(H::Matrix)

  el = eltype(H[1])
  one = el(1)

  tol = 1000*eps(abs.(one))

  r,c = size(H)
  if r!=c
    display("H needs to be a square matrix: size(H)= $r,$c")
  end

  W  = Matrix{vt}(I,r,c)
  V  = Matrix{vt}(I,r,c)
  I0 = Matrix{vt}(I,r,c)

  one = vt(1)

  for i in 1:c-2

    a = vt(1)
    b = vt(0)
    d = vt(1)
    c = -H[i+2,i]/H[i+1,i]
    if abs(H[i+1,i])<tol
      x3a = abs(H[i+1,i])
      println("Possible Division by 0 abs(x3)=$x3a")
    end  
   
#   Left Multiplication      
    Q       = copy(I0)
    Q[i+2,i+1] = c
    A       = Q*H
    W       = Q*W

#   Right multiplication
    Q[i+2,i+1] = -c
    H      .= A*Q
    V       = V*Q
    
  end

  return V,W
end
#----------------------------------------------------------------------
function ChaseBulgeTriDiagonal2!(H::Matrix)

# Modify individual entries
# I assume the the bulge size is 1

  el = eltype(H[1])
  one = el(1)

  tol = 1000*eps(abs.(one))

  row,col = size(H)
  if row!=col
    display("H needs to be a square matrix: size(H)= $r,$c")
  end

  W  = Matrix{vt}(I,row,col)
  V  = Matrix{vt}(I,row,col)
  I0 = Matrix{vt}(I,row,col)

  one = vt(1)
  zro = vt(0)

  for i in 1:col-2

    x1 = H[i+1,i]
    x2 = H[i+2,i]
    y1 = H[i,i+1]
    y2 = H[i+2,i+1]
    z1 = H[i+1,i+2]
    if (i<col-2)
      z2 = H[i+3,i+2]
    else
      z2 = zro
    end
    if i+3<=col
      w1 = H[i+2,i+3]
    else
      w1 = zro
    end  
   
    α    = H[i,i]
    β    = H[i+1,i+1]
    γ    = H[i+2,i+2]

#   Ensure that the determinant is one.
#   Otherwise we need to factor that in the inverse
    a = one
    b = zro
    d = one
    c = -d*x2/x1
    if abs(x1)<tol
      x1a = abs(x1)
      println("Possible Division by 0 abs(x1)=$x1a")
    end  
   
#   Left Multiplication      
    Q       = copy(I0)
    Q[i+2,i+1] = c
    
    W       = Q*W

#   Right multiplication
    Q[i+2,i+1] = -c
    V          = V*Q

    H[i+1,i]      = a*x1
    H[i+2,i]      = zro

    H[i,i+1]      = d*y1
    H[i+1,i+1]    = a*β*d - a*c*z1
    H[i+2,i+1]    = d*(c*β + d*y2) - c*(c*z1 + γ*d)

    H[i+2,i+2]    = a*(c*z1 + γ*d)
    if (i<col-2)
      H[i+3,i+1]  = -c*z2
      H[i+3,i+2]  = a*z2
    end  

    if i+3<=col
      H[i+2,i+3]  = d*w1
    end  
   
  end

  return V,W
end
#----------------------------------------------------------------------
function ChaseBulgeTriDiagonal3!(H::Matrix)

# Chase Upper Bulge  
# Modify individual entries
# I assume the the bulge size is 1

  el = eltype(H[1])
  one = el(1)

  tol = 1000*eps(abs.(one))

  row,col = size(H)
  if row!=col
    display("H needs to be a square matrix: size(H)= $r,$c")
  end

  W  = Matrix{vt}(I,row,col)
  V  = Matrix{vt}(I,row,col)
  I0 = Matrix{vt}(I,row,col)

  one = vt(1)
  zro = vt(0)

  for i in 1:col-2

    x1 = H[i,i+1]
    x2 = H[i,i+2]  

    y1 = H[i+1,i]
    y2 = H[i+1,i+2]

    z1 = H[i+2,i+1]
    if (i+2<col)
      z2 = H[i+2,i+3]
    else
      z2 = zro
    end
    if i+3<=col
      w1 = H[i+3,i+2]
    else
      w1 = zro
    end  

    α    = H[i,i]
    β    = H[i+1,i+1]
    γ    = H[i+2,i+2]

#   Ensure that the determinant is one.
#   Otherwise we need to factor that in the inverse
    a = one
    c = zro
    d = one
    b = a*x2/x1
    if abs(x1)<tol
      x1a = abs(x1)
      println("Possible Division by 0 abs(x1)=$x1a")
    end  
   
#   Left Multiplication      
    Q          = copy(I0)
    Q[i+1,i+1] = a
    Q[i+1,i+2] = b
    Q[i+1,i+1] = d
   
    W       = Q*W

#   Right multiplication
    Q[i+1,i+1] = d
    Q[i+1,i+2] = -b
    Q[i+1,i+1] = a

    V       = V*Q

    H[i,i+1]      = d*x1
    H[i,i+2]      = zro
   
    H[i+1,i]      = a*y1
    H[i+1,i+1]    = a*β*d + d*b*z1
    H[i+1,i+2]    = a*(a*y2 + b*γ) - b*(b*z1 + a*β)
    if (i+2<col)
      H[i+1,i+3]  = b*z2
      H[i+2,i+3]  = d*z2
    end

    H[i+2,i+1]    = d*d*z1
    H[i+2,i+2]    = (a*γ - b*z1)*d

    if i+3<=col
      H[i+3,i+2]  = a*w1
    end  
    
  end

  return V,W
end
#----------------------------------------------------------------------











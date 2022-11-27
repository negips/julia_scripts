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
  el    = eltype(T[1])
  one   = el(1)
  tol   = 1000.0*eps(abs(one))
  αtol   = 1.0e-2

  V     = Matrix{el}(I,r,c)
  W     = Matrix{el}(I,r,c)

## Create a bulge in the matrix
#  if abs(T[2,1])<tol
#    Wi,Vi = SmallX1_fix!(T,1)
##   Collect Right multipliers 
#    V = V*Vi
##   Collect Left multipliers
#    W = Wi*W
#  end  

#  T,v,w     = CreateLowerBulgeOblique(T,λ)
  if  (abs(T[1,1]-λ) > αtol && abs(T[2,2]-λ) > αtol)
    if (abs(T[2,1])< tol)
      Vi,Wi = SmallX1_fix!(T)
      V = V*Vi
      W = Wi*W
    else  
      w,v       = SimilarityTransformBulge!(T,λ,1)
    end  
    Vi,Wi     = ChaseBulgeTriDiagonal2!(T,λ)
  else
    println("Small α. Hybrid Chase")
    T,Q       = CreateBulge(T,1,λ,1)    
    Vi,Wi     = HybridBulgeChase!(T)  
  end  

# Collect Right multipliers 
  V = V*Vi
# Collect Left multipliers
  W = Wi*W

  return T,V,W
end
#----------------------------------------------------------------------
# Lower Bulge Chase Algorithm with Oblique projectors
function NegiAlgHybrid(T::Matrix,λ)
# Also known as BulgeChase Algorithm  

  # H       - Tridiagonal Matrix
  # λ       - Shift

  r,c   = size(T)
  el    = eltype(T[1])
  one   = el(1)
  tol   = 1000.0*eps(abs(one))
  αtol   = 1.0e-2

  V     = Matrix{el}(I,r,c)
  W     = Matrix{el}(I,r,c)

## Create a bulge in the matrix
#  if abs(T[2,1])<tol
#    Wi,Vi = SmallX1_fix!(T,1)
##   Collect Right multipliers 
#    V = V*Vi
##   Collect Left multipliers
#    W = Wi*W
#  end  

#  T,v,w     = CreateLowerBulgeOblique(T,λ)
  T,Q       = CreateBulge(T,1,λ,1)    
  Vi,Wi     = HybridBulgeChase!(T)  

# Collect Right multipliers 
  V = V*Vi
# Collect Left multipliers
  W = Wi*W

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
# Lower Right Bulge Chase (up) Algorithm with Oblique projectors
function NegiAlg4(T::Matrix,λ)
# Also known as BulgeChase Algorithm  

  # H       - Tridiagonal Matrix
  # λ       - Shift

  r,c   = size(T)

# Create a bulge in the matrix
  T,v,w     = CreateLowerRightBulgeOblique(T,λ)
  Vi,Wi     = ChaseBulgeTriDiagonal4!(T)

# Collect Left multipliers 
  V = Vi
# Collect Right multipliers
  W = Wi

  return T,V,W
end
#----------------------------------------------------------------------

function CreateLowerBulgeOblique(T::Matrix,λ)

# Create bulge in the sub diagonal  

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
  T        .= T + λ*I

  return T,x,y
end

#----------------------------------------------------------------------
function CreateLowerBulgeOblique2(T::Matrix,λ,k::Int)

# For creating a new bulge in the middle of the 
# matrix

  # T       - Tri-diagonal Matrix
  # λ       - Shifts

  T         = T - λ*I
  el        = eltype(T[1])
  zro       = el(0)
  one       = el(1)
  r,c       = size(T)
  x         = zeros(el,r)
  x[k]      = zro
  x[k+1]    = T[k+1,k]
  y         = transpose(zeros(vt,r)) 
  y[k:c]    = transpose(T[k,k:c])
  α         = (y*x)[1]
  y         = y/α      # This makes <y,x> = 1.0
  β         = (y*T[:,k])[1]
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

# For creating a bulge in the super diagonal elements  

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
function CreateLowerRightBulgeOblique(T::Matrix,λ)

# For creating a bulge in the subdiagonal 
# on the lower right side of the matrix

  # T       - Tri-diagonal Matrix
  # λ       - Shifts

  T         = T - λ*I
  el        = eltype(T[1])
  zro       = el(0)
  one       = el(1)
  r,c       = size(T)
  y         = transpose(zeros(el,c))
  y[c-1]    = T[r,c-1]

  x         = T[:,c]
 
  α         = (y*x)[1]
  x         = x/α      # This makes <y,x> = 1.0
  β         = (transpose(T[r,:])*x)[1]
  x         = x/β      # => <y,x> = 1.0/β
  V         = I - (x*y)
  Vi        = I + (x*y)/(one - one/β)

# This creates the Bulge  
# A = Q0*H*Q0  
  tmp       = T*V
  T        .= Vi*tmp
  T         = T + λ*I

  return T,x,y
end

#----------------------------------------------------------------------
function CreateUpperRightBulgeOblique(T::Matrix,λ)

# Unfortunate nomenclature.
# Creates a bulge in the 2nd superdiagonal
# at the bottom of the matrix

  # T       - Tri-diagonal Matrix
  # λ       - Shifts

  T         = T - λ*I
  el        = eltype(T[1])
  zro       = el(0)
  one       = el(1)
  r,c       = size(T)

  tol       = 1000*eps(abs(one))

  V          = Matrix{vt}(I,r,c)
  Vi         = Matrix{vt}(I,r,c)
  cs         = T[r,c]
  sn         = -T[r-1,c]
  rad        = sqrt(cs*cs + sn*sn)
  V[r-1,c-1] = cs
  V[r-1,c]   = -sn
  V[r,c]     = 1.0/cs

  Vi[r-1,c-1] = one/cs
  Vi[r-1,c]   = sn
  Vi[r,c]     = cs
    
  x           = zeros(vt,r)
  x[r]        = cs
  x[r-1]      = -sn    

  y           = transpose(zeros(vt,c))
  y[c]        = one

# This creates the Bulge  
# A = Q0*H*Q0  
  tmp       = T*V
  T        .= Vi*tmp
  T         = T + λ*I

  return T,x,y
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
function HybridBulgeChase!(H::Matrix)

  el = eltype(H[1])
  one = el(1)

  tol = 100*eps(abs.(one))

  r,c = size(H)
  if r!=c
    display("H needs to be a square matrix: size(H)= $r,$c")
  end

  W  = Matrix{vt}(I,r,c)
  V  = Matrix{vt}(I,r,c)
  I0 = Matrix{vt}(I,r,c)

  for i in 1:c-2

    Q       = ChaseBulgeDownOneStep!(H,i)  
#   Left Multiplication      
    W       = Q*W
#   Right multiplication
    V       = V*Q

    Vi,Wi   = ChaseUpperBulgeObliqueOneStep!(H,i)
    V       = V*Vi
    W       = Wi*W
    
  end

  return V,W
end
#----------------------------------------------------------------------

function ChaseBulgeTriDiagonal!(H::Matrix)

  el = eltype(H[1])
  one = el(1)

  tol = 100*eps(abs.(one))

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
function ChaseBulgeTriDiagonal2!(H::Matrix,λ)

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

    if abs(x1)<1.0e-4  #tol
      x1a = abs(x1)
      x2a = abs(x2)
#      println("Possible Division by 0 in Chase2 abs(x1)=$x1a, abs(x2)=$x2a, i=$i")
#      Ql,Qr = SmallX1_fix!(H,i)
#      Ql,Qr = SimilarityTransformBulge!(H,λ,i+1)     
#      V = V*Qr
#      W = Ql*W
    end

    if abs(x2)<tol
#      println("Recreating Reflector")
#      Ql,Qr = SmallX1_fix!(H,i)
#      V = V*Qr
#      W = Ql*W
    end
      
    Ql,Qr   = SimilarityTransform!(H,i,col)
    V = V*Qr
    W = Ql*W
   
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
      println("Possible Division by 0 in Chase3 abs(x1)=$x1a")
    end  
   
#   Left Multiplication      
    Q          = copy(I0)
    Q[i+1,i+1] = a
    Q[i+1,i+2] = b
    Q[i+2,i+2] = d
   
    W       = Q*W

#   Right multiplication
    Q[i+1,i+1] = d
    Q[i+1,i+2] = -b
    Q[i+2,i+2] = a

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
function ChaseBulgeTriDiagonal4!(H::Matrix)

# Chase Bottom Left to top 
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

  for i in col:-1:3

    x1 = H[i,i-1]
    x2 = H[i,i-2]  

    y1 = H[i-1,i]
    y2 = H[i-1,i-2]

    z1 = H[i-2,i-1]
    if (i-2>1)
      z2 = H[i-2,i-3]
    else
      z2 = zro
    end

    if i-3>1
      w1 = H[i-3,i-2]
    else
      w1 = zro
    end  

    α    = H[i,i]
    β    = H[i-1,i-1]
    γ    = H[i-2,i-2]

#   Ensure that the determinant is one.
#   Otherwise we need to factor that in the inverse
    a = one
    b = zro
    d = one
    c = d*x2/x1
    if abs(x1)<tol
      x1a = abs(x1)
      println("Possible Division by 0 in Chase4, abs(x1)=$x1a")
    end  
   
#   Left Multiplication      
    Q          = copy(I0)
    Q[i-1,i-1] = d
    Q[i-1,i-2] = c
    Q[i-2,i-2] = a
   
    W       = Q*W

#   Right multiplication
    Q[i-1,i-1] = a
    Q[i-1,i-2] = -c
    Q[i-2,i-2] = d


    V       = V*Q

    H[i,i-1]      = a*x1
    H[i,i-2]      = zro
   
    H[i-1,i]      = d*y1
    H[i-1,i-1]    = a*β*d + a*c*z1
    H[i-1,i-2]    = d*(d*y2 + c*γ) - c*(c*z1 + d*β)
    if (i-2>1)
      H[i-1,i-3]  = c*z2
      H[i-2,i-3]  = a*z2
    end

    H[i-2,i-1]    = a*a*z1
    H[i-2,i-2]    = a*(d*γ - c*z1)

    if i-3>1
      H[i-3,i-2]  = d*w1
    end  
    
  end

  return V,W
end
#----------------------------------------------------------------------

function ChaseBulgeTriDiagonal5!(H::Matrix)

# Chase Bottom Right to top 
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

  for i in col:-1:3

    x1 = H[i-1,i]
    x2 = H[i-1,i]  

    y1 = H[i,i-1]
    y2 = H[i-2,i-1]

    z1 = H[i-1,i-2]
    if (i-2>1)
      z2 = H[i-3,i-2]
    else
      z2 = zro
    end

    if i-3>1
      w1 = H[i-2,i-3]
    else
      w1 = zro
    end  

    α    = H[i,i]
    β    = H[i-1,i-1]
    γ    = H[i-2,i-2]

#   Ensure that the determinant is one.
#   Otherwise we need to factor that in the inverse
    a = one
    b = zro
    d = one
    c = -d*x2/x1
    if abs(x1)<tol
      x1a = abs(x1)
      println("Possible Division by 0 in Chase5, abs(x1)=$x1a")
    end  
   
#   Left Multiplication      
    Q          = copy(I0)
    Q[i-1,i-1] = d
    Q[i-1,i-2] = c
    Q[i-2,i-2] = a
   
    W       = Q*W

#   Right multiplication
    Q[i-1,i-1] = a
    Q[i-1,i-2] = -c
    Q[i-2,i-2] = d


    V       = V*Q

    H[i-1,i]      = a*x1
    H[i-2,i]      = zro
   
    H[i,i-1]      = d*y1
    H[i-1,i-1]    = a*β*d - a*c*z1
    H[i-2,i-1]    = d*(d*y2 + c*β) - c*(c*z1 + d*γ)
    if (i-2>1)
      H[i-3,i-1]  = -c*z2
      H[i-3,i-2]  = a*z2
    end

    H[i-1,i-2]    = a*a*z1
    H[i-2,i-2]    = a*(d*γ + c*z1)

    if i-3>1
      H[i-2,i-3]  = d*w1
    end  
    
  end

  return V,W
end
#----------------------------------------------------------------------
function SimilarityTransform!(H::Matrix,i::Int,n::Int)

#   ( I  0  0  0 )   ( α   y1  0   0  )   (I   0   0   0 )
#   ( 0  a  b  0 ) X ( x1  β   z1  0  ) X (0   a  -b   0 )
#   ( 0  c  d  0 )   ( x2  y2  γ   w1 )   (0  -c   d   0 ) X 1/(ad - bc)
#   ( 0  0  0  I )   ( 0   0   z2  δ  )   (0   0   0   I )

    el  = eltype(H[1])
    zro = el(0)
    one = el(1)

    tol = 1000*eps(abs(one))
    I0 = Matrix{el}(I,n,n)

#   Left Multiplication
    α    = H[i,i]
    β    = H[i+1,i+1]
    γ    = H[i+2,i+2]
  
    x1 = H[i+1,i]
    x2 = H[i+2,i]
    y1 = H[i,i+1]
    y2 = H[i+2,i+1]
    z1 = H[i+1,i+2]
    if (i+3<=n)
      z2 = H[i+3,i+2]
      δ  = H[i+3,i+3]
      w1 = H[i+2,i+3]
      y3 = H[i+3,i+1]
    else
      z2 = zro
      δ  = zro
      w1 = zro
      y3 = zro
    end

    if abs(x1)<tol
      x1a = abs(x1)
      x2a = abs(x2)
#      println("Division by 0 abs(x1)=$x1a, abs(x2)=$x2a, i=$i")
    end

    b = zro
    d = one
    c = -d*x2/x1
    a = one
    D0 = a*d - b*c
  
    H[i,i]        = α
    H[i+1,i]      = a*x1 + b*x2
#    H[i+2,i]      = c*x1 + d*x2
    H[i+2,i]      = zro       # Forcing exact zero      

    H[i,i+1]      = y1
    H[i+1,i+1]    = a*β + b*y2
    H[i+2,i+1]    = c*β + d*y2

    H[i+1,i+2]    = a*z1 + b*γ
    H[i+2,i+2]    = c*z1 + γ*d

    if (i+3<=n)
      H[i+1,i+3]  = b*w1
      H[i+2,i+3]  = d*w1
    end  

    Ql            = copy(I0)
    Ql[i+1,i+1]   = a
    Ql[i+2,i+1]   = c
    Ql[i+1,i+2]   = b
    Ql[i+2,i+2]   = d
#----------------------------------------

#   Right Multiplication
    α    = H[i,i]
    β    = H[i+1,i+1]
    γ    = H[i+2,i+2]
  
    x1   = H[i+1,i]
    x2   = H[i+2,i]
    y1   = H[i,i+1]
    y2   = H[i+2,i+1]
    z1   = H[i+1,i+2]
    if (i+3<=n)
      z2 = H[i+3,i+2]
      δ  = H[i+3,i+3]
      w1 = H[i+2,i+3]
      y3 = H[i+3,i+1]
    else
      z2 = zro
      δ  = zro
      w1 = zro
    end

    dold = d
    b    = -b/D0
    c    = -c/D0
    d    = a/D0
    a    = dold/D0

#   Right multiplication
    H[i,i+1]      = a*y1
    H[i+1,i+1]    = a*β  + c*z1
    H[i+2,i+1]    = a*y2 + c*γ

    H[i,i+2]      = b*y1
    H[i+1,i+2]    = d*z1 + b*β
    H[i+2,i+2]    = d*γ  + b*y2

    if (i+3<=n)
      H[i+3,i+1]  = a*y3 + c*z2
      H[i+3,i+2]  = d*z2 + b*y3
    end  

    Qr            = copy(I0)
    Qr[i+1,i+1]   = a
    Qr[i+2,i+1]   = c
    Qr[i+1,i+2]   = b
    Qr[i+2,i+2]   = d

   return Ql,Qr
end
#---------------------------------------------------------------------- 
function SmallX1_fix!(H::Matrix)

#   ( I  0  0  0 )   ( α   y1  0   0  )   ( I    0    0   0 )
#   ( b  a  0  0 ) X ( x1  β   z1  0  ) X ( -b   1    0   0 )
#   ( 0  0  I  0 )   ( x2  y2  γ   w1 )   ( 0    0    I   0 ) X 1/a
#   ( 0  0  0  I )   ( 0   0   z2  δ  )   ( 0    0    0   I )

    i = 1

    el  = eltype(H[1])
    n,n1 = size(H)
    zro = el(0)
    one = el(1)

    I0 = Matrix{el}(I,n,n)

#   Left Multiplication
    α    = H[i,i]
    β    = H[i+1,i+1]
    γ    = H[i+2,i+2]
  
    x1 = H[i+1,i]
    x2 = H[i+2,i]
    y1 = H[i,i+1]
    y2 = H[i+2,i+1]
    z1 = H[i+1,i+2]

    b = one
    a = one       # Diagonal element. Can't be zero
    D0 = a
  
    H[i,i]        = α
    H[i+1,i]      = a*x1 + b*α

    H[i+1,i+1]    = a*β + b*y1

    H[i+1,i+2]    = a*z1

    Ql            = copy(I0)
    Ql[i+1,i+1]   = a
    Ql[i+1,i]     = b
#----------------------------------------

#   Right Multiplication
    α    = H[i,i]
    β    = H[i+1,i+1]
    γ    = H[i+2,i+2]
  
    x1 = H[i+1,i]
    x2 = H[i+2,i]
    y1 = H[i,i+1]
    y2 = H[i+2,i+1]
    z1 = H[i+1,i+2]

    b = one
    a = one       # Diagonal element. Can't be zero
    D0 = a
  
    H[i,i]        = α
    H[i+1,i]      = a*x1 + b*α

    H[i+1,i+1]    = a*β + b*y1

    H[i+1,i+2]    = a*z1

    b    = -b/D0
    a    =  1/D0

#   Right multiplication
    H[i,i]        = α + b*y1     
    H[i+1,i]      = b*β + x1
    H[i+2,i]      = b*y2 + x2

    H[i,i+1]      = a*y1
    H[i+1,i+1]    = a*β
    H[i+2,i+1]    = a*y2

    Qr            = copy(I0)
    Qr[i+1,i+1]   = a
    Qr[i+1,i]     = b

   return Ql,Qr
end
#----------------------------------------------------------------------
function CreateBulgeMiddle!(H::Matrix,λ,i::Int)

#   ( I  0  0  0 )   ( δ   x1  0   0  )   ( I    0    0   0 )
#   ( 0  I  0  0 ) X ( w2  α   y1  0  ) X ( 0    1    0   0 )
#   ( 0  c  I  0 )   ( 0   x2  β   z1 )   ( 0   -c    I   0 )
#   ( 0  0  0  I )   ( 0   0   y2  γ  )   ( 0    0    0   I )

#   Use this when x2 (above) is small

#   The index "i" here corresponds to the H[i,i] == α

    el  = eltype(H[1])
    n,n1 = size(H)
    zro = el(0)
    one = el(1)

    I0 = Matrix{el}(I,n,n)

#   Left Multiplication
    α    = H[i,i]       - λ
    β    = H[i+1,i+1]   - λ
    γ    = H[i+2,i+2]   - λ

    if (i>1)
      x1 = H[i-1,i]
      w2 = H[i,i-1]
    else
      x1 = zro
      w2 = zro
    end  

    x2 = H[i+1,i]
    y1 = H[i,i+1]
    y2 = H[i+2,i+1]
    if i+2<=n 
      z1 = H[i+1,i+2]
    else
      z1 = zro
    end  

    a  = one
    c  = one
    d  = one
    D0 = a*d

    if (i>1)
      H[i,i-1]    = a*w2
      H[i+1,i-1]  = c*w2
    end  

    H[i,i]        = a*α
    H[i+1,i]      = c*α  + d*x2

    H[i,i+1]      =        a*y1
    H[i+1,i+1]    = d*β  + c*y1

    if i+2<=n
      H[i+1,i+2]  =        d*z1
      H[i+2,i+2]  = γ
    end

    Ql            = copy(I0)
    Ql[i+1,i+1]   = a
    Ql[i+1,i+1]   = d
    Ql[i+1,i]     = c
#----------------------------------------
#   Right Multiplication
    α    = H[i,i]    
    β    = H[i+1,i+1]
    γ    = H[i+2,i+2]

    if (i>1)
      x1 = H[i-1,i]
      w2 = H[i,i-1]
      w3 = H[i+1,i-1]
    else
      x1 = zro
      w2 = zro
      w3 = zro
    end  

    x2 = H[i+1,i]
    y1 = H[i,i+1]
    y2 = H[i+2,i+1]
    if i+2<=n 
      z1 = H[i+1,i+2]
    else
      z1 = zro
    end  

    aold = a
    a  = d/D0
    c  = -c/D0
    d  = aold/D0

    if (i>1)
      H[i-1,i]    = a*x1
    end  

    H[i,i]        = a*α  + c*y1
    H[i+1,i]      = c*β  + a*x2

    H[i,i+1]      =        d*y1
    H[i+1,i+1]    = d*β

    if i+2<=n
      H[i+1,i+2]  =          z1
      H[i+2,i+2]  = γ
     
      H[i+2,i]    =        c*y2
      H[i+2,i+1]  =        d*y2
    end  

#   Remove shifts
    H[i,i]        = H[i,i]     + λ  
    H[i+1,i+1]    = H[i+1,i+1] + λ  
    H[i+2,i+2]    = H[i+2,i+2] + λ  

    Qr            = copy(I0)
    Qr[i+1,i+1]   = a
    Qr[i+1,i+1]   = d
    Qr[i+1,i]     = c

   return Ql,Qr
end
#---------------------------------------------------------------------- 

function SimilarityTransformBulge!(H::Matrix,λ,i::Int)

#   Create bulge in the sub diagonal via similarity transform
#   ( a  0  0 )   ( α   y1  0  )   (d   0  0 )
#   ( c  d  0 ) X ( x1  β   z1 ) X (-c  a  0 ) X 1/(ad - 0)
#   ( 0  0  I )   ( 0   y2  γ  )   (0   0  I )

    el  = eltype(H[1])
    zro = el(0)
    one = el(1)
    tol = 1.0e-2

    I0 = Matrix{el}(I,size(H))

#   Left Multiplication
    α    = H[i,i]       - λ
    β    = H[i+1,i+1]   - λ
    γ    = H[i+2,i+2]   - λ
  
    x1 = H[i+1,i]
    x2 = H[i+2,i]       # This should be zero
    y1 = H[i,i+1]
    y2 = H[i+2,i+1]
    z1 = H[i+1,i+2]

    if (i>1)
      w2 = H[i,i-1]
    end  

    b = zro
    d = one
    c = -d*x1/α
    a = one
    D0 = a*d - b*c
  
    H[i,i]        = a*α
    H[i+1,i]      = zro       # Forcing exact zero

    H[i,i+1]      = a*y1
    H[i+1,i+1]    = c*y1 + d*β

    H[i+1,i+2]    = d*z1

    H[i+2,i+2]    = γ

    if (i>1)
      H[i,i-1]    = a*w2
      H[i+1,i-1]  = c*w2
    end  

    Ql            = copy(I0)
    Ql[i,i]       = a
    Ql[i+1,i]     = c
    Ql[i,i+1]     = b
    Ql[i+1,i+1]   = d
#----------------------------------------

#   Right Multiplication
    α    = H[i,i]
    β    = H[i+1,i+1]
    γ    = H[i+2,i+2]
  
    x1 = H[i+1,i]
    x2 = H[i+2,i]       # This should be zero
    y1 = H[i,i+1]
    y2 = H[i+2,i+1]
    z1 = H[i+1,i+2]

    dold = d
    b    = zro
    c    = -c/D0
    d    = a/D0
    a    = dold/D0

#   Right multiplication
    H[i,i]        = a*α  + c*y1
    H[i+1,i]      = c*β
    H[i+2,i]      = c*y2

    H[i,i+1]      = d*y1
    H[i+1,i+1]    = d*β
    H[i+2,i+1]    = d*y2

    H[i+2,i+2]    = γ

#   Remove shifts
    H[i,i]        = H[i,i]     + λ  
    H[i+1,i+1]    = H[i+1,i+1] + λ  
    H[i+2,i+2]    = H[i+2,i+2] + λ  

    Qr            = copy(I0)
    Qr[i,i]       = a
    Qr[i+1,i]     = c
    Qr[i,i+1]     = b
    Qr[i+1,i+1]   = d

    return Ql,Qr

end

#----------------------------------------------------------------------





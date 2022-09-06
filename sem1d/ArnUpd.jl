# Arnoldi Update
function ArnUpd(V::Matrix,b::Int,B::Vector,v::Vector,k::Int,ngs::Int)

# V         - Krylov Vector
# H         - Upper Hessenberg
# b         - Block Size
# B         - Weight (vector)
# β         - previous/new residual norm
# k         - Current Krylov size
# kmax      - Maximum Krylov size
# v         - New vector v = Ax
# ngs       - No of Gram-Schmidt


  h   = zeros(eltype(v),k)
  a   = real(v[1])
  localprec = typeof(a)
  β   = localprec(0)
  fac = localprec(1000)

#  if (typeof(a)==Float64)
#    β   = 0.
#    fac = 1000.0
#  else
#    β   = BigFloat(0.)
#    fac = BigFloat(1000.0)
#  end  
 
  r = copy(v)

  tol = eps(fac)

# New Arnoldi Vector
  if k > 0
    h = V[:,1:k]'*(B.*r)
    q  = similar(r)
    mul!(q,V[:,1:k],h)
    r  .= r .- q
    if ngs>1
      for j in 2:ngs
        g = V[:,1:k]'*(B.*r)
        mul!(q,V[:,1:k],g)
        r .= r .- q
        h .= h .+ g
      end
    elseif ngs <=0 
      while true
        g   = V[:,1:k]'*(B.*r)
        res = abs(g'*g)
        if res<tol
          break
        end
        h = h .+ g
        r = r - V[:,1:k]*g
      end  
    end
    β        = norm(sqrt(r'*(B.*r)))
    r        = r/β
  else
    β        = norm(sqrt(r'*(B.*r)))
    r        = r/β
  end

  return h,β,r

end

#----------------------------------------------------------------------
function ArnUpd(V::Matrix,b::Int,B::Matrix,v::Vector,k::Int,ngs::Int)

# V         - Krylov Vector
# H         - Upper Hessenberg
# b         - Block Size
# B         - Weight (Matrix)
# β         - previous/new residual norm
# k         - Current Krylov size
# kmax      - Maximum Krylov size
# v         - New vector v = Ax
# ngs       - No of Gram-Schmidt


  h = zeros(eltype(v),k)
  a   = real(v[1])
  localprec = typeof(a)
  β   = localprec(0)
  fac = localprec(1000)

#  if (typeof(a)==Float64)
#    β   = 0.
#    fac = 1000.0
#  else
#    β   = BigFloat(0.)
#    fac = BigFloat(1000.0)
#  end  
  r = copy(v)

  tol = eps(fac)

# New Arnoldi Vector
  if k >= b
    h = V[:,1:k]'*(B*r)
    q  = similar(r)
    mul!(q,V[:,1:k],h)
    r  .= r .- q
    if ngs>1
      for j in 2:ngs
        g = V[:,1:k]'*(B*r)
        mul!(q,V[:,1:k],g)
        r .= r .- q
        h .= h .+ g
      end
    elseif ngs <=0 
      while true
        g   = V[:,1:k]'*(B*r)
        res = abs(g'*g)
        if res<tol
          break
        end
        h = h .+ g
        r = r - V[:,1:k]*g
      end  
    end
    β        = norm(sqrt(r'*(B*r)))
    r        = r/β
  else
    β        = norm(sqrt(r'*(B*r)))
    r        = r/β
  end

  return h,β,r

end   
#----------------------------------------------------------------------



# Arnoldi Update
function ArnUpd(V::Matrix,b::Vector,v::Vector,k::Int,ngs::Int)

# V         - Krylov Vector
# H         - Upper Hessenberg
# b         - Weight (vector)
# β         - previous/new residual norm
# k         - Current Krylov size
# kmax      - Maximum Krylov size
# v         - New vector v = Ax
# ngs       - No of Gram-Schmidt


  h = zeros(typeof(v[1]),k)
  β   = 0.
  r = copy(v)

# New Arnoldi Vector
  if k > 0
    h = V[:,1:k]'*(b.*r)
    q  = similar(r)
    mul!(q,V[:,1:k],h)
    r  .= r .- q
    if ngs>1
      for j in 2:ngs
        g = V[:,1:k]'*(b.*r)
        mul!(q,V[:,1:k],g)
        r .= r .- q
        h .= h .+ g
      end  
    end
    β        = norm(sqrt(r'*(b.*r)))
    r        = r/β
  else
    β           = norm(sqrt(r'*(b.*r)))
    r           = r/β
  end

  return h,β,r

end   

# Arnoldi Update
function ArnUpd!(V::Matrix,H::Matrix,b::Vector,k::Int,kmax::Int,v::Vector,ngs::Int)
  global V,H,v

# Update Arnoldi Vector
  if k > 0
    h = V[:,k]'*(b.*v)
    v = v .- V[:,1:k]*h
    if ngs>1
      for j in 2:ngs
        g = V[:,1:k]'*(b.*v)
        v = v .- V[:,1:k]*g
        h = h .+ g
      end  
    end
    H[1:k,k] = h
    β                = sqrt(v'*(b.*v))
    H[nkryl+1,nkryl] = β
    v                = v/β
    k                = k + 1
    if k<=kmax
      V[:,k]       = v
    end
  else
    β           = sqrt(v'*(b.*v))
    v           = v/β
    V[:,1]      = v
    k           = 1
  end  

  return V,H,v,k
end   

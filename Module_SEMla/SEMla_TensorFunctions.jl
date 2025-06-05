#     Non-Constructor Tensor-Function Definitions
#---------------------------------------------------------------------- 
@doc raw"""

      function tensorfieldOP!(tfout::AbstractArray{T,N},tfin::AbstractArray{T,N},OP::AbstractMatrix{T},wk::AbstractArray{T,N}) where {T<:Number}
        
      Operation (OP) on two-dimensional tensor arrays (tfin).
      OP is applied along all dimensions.

"""
@views function tensorfieldOP!(to::AbstractTensorField2{T},ti::AbstractTensorField2{T},OP::AbstractMatrix{T},wk::AbstractMatrix{T}) where {T<:Number}

  nel      = ti.nel
  p1,p2    = ti.n
  q1,q2    = to.n
  r1,r2    = size(OP)

  OPT      = OP'

  for e in 1:nel
    # These should be views            
    w1     = view(wk, 1:r1,1:p2)
    v1     = view(ti.tfield, 1:p1,1:p2,e)
    u1     = view(to.tfield, 1:r1,1:r1,e)

    # wk = OP*v1
    mul!(w1,OP,v1)

    # to = wk*(OP^T)
    mul!(u1,w1,OPT)
  end  

  return nothing
end
#----------------------------------------------------------------------       
@doc raw"""

      function tensorfieldOP12!(tfout::AbstractArray{T,N},tfin::AbstractArray{T,N},OP1::AbstractMatrix{T},OP2::AbstractMatrix{T},wk::AbstractArray{T,N}) where {T<:Number}
        
      Operation (OP1 and OP2) on two-dimensional tensor arrays (tfin).
      OP1 along the first dimension and OP2 along the second dimension.

"""
@views function tensorfieldOP12!(to::AbstractTensorField2{T},ti::AbstractTensorField2{T},OP1::AbstractMatrix{T},OP2::AbstractMatrix{T},wk::AbstractMatrix{T}) where {T<:Number}


  nel      = ti.nel
  p1,p2    = ti.n
  q1,q2    = to.n
  r1,r2    = size(OP1)
  s1,s2    = size(OP2)

  OP2T     = OP2'

  for e in 1:nel
    # These should be views            
    w1       = view(wk, 1:r1,1:p2)
    v1       = view(ti.tfield, 1:p1,1:p2,e)
    u1       = view(to.tfield, 1:r1,1:s1,e)

    # wk = OP1*v1
    mul!(w1,OP1,v1)

    # u1 = wk*(OP2^T)
    mul!(u1,w1,OP2T)
  end  

  return nothing
end
#----------------------------------------------------------------------       
@doc raw"""

      function tensorOPn!(tfout::AbstractArray{T},tfin::AbstractArray{T},OP::AbstractMatrix{T},n::Int) where {T<:Number}
          
      Operation on the tensor arrays along dimension n.
      If n>1 it is assumed that an adjoint operator has been passed

"""
@views function tensorOPn!(tfout::AbstractArray{T,N},tfin::AbstractArray{T,N},OP::AbstractMatrix{T},n::Int) where {T<:Number,N}

# n is the dimension along which the tensor operator OP acts.

  si       = size(tfin)
  so       = size(tfout)
  sop      = size(OP)
  
  if n == 1
    # We should have the Direct operator 
    @assert si[1] == sop[2]    "OP and Input Tensor dimension mismatch $sop, $si"
    i1  = si[n]
    j1  = prod(si[n+1:end])
    fi  = reshape(tfin,(i1,j1))
  
    i2  = so[n]
    j2  = prod(so[n+1:end])
    fo  = reshape(tfout,(i2,j2))
  
    ci = CartesianIndices((1:i1,1:j1))
    co = CartesianIndices((1:i2,1:j1))
  
    mul!(fo[co],OP,fi[ci])
   
  else
    # We should have the Adjoint operator 
    @assert si[n] == sop[1]    "OP and Input Tensor dimension mismatch $sop, $si"
    i1  = prod(si[1:n-1])
    j1  = si[n]
    k1  = prod(si[n+1:end])
    fi  = reshape(tfin,(i1,j1,k1))
  
    i2  = prod(so[1:n-1])
    j2  = so[n]
    k2  = prod(so[n+1:end])
    fo  = reshape(tfout,(i2,j2,k2))
  
    ci = CartesianIndices((1:i1,1:j1))
    co = CartesianIndices((1:i2,1:j2))
   
    for kk in 1:k1
      mul!(fo[co,kk],fi[ci,kk],OP)
    end 
  
  end
  
  return nothing
end
#----------------------------------------------------------------------       
@doc raw"""

      tensordimension(AbstractTensorField{T,N}) where{T:<Number,N}

      Returns the tensor field dimension: N-1

"""
function tensordimension(f::AbstractTensorField{T,N}) where{T<:Number,N}
  return N-1
end

#----------------------------------------------------------------------       
@doc raw"""

      tensortype(AbstractTensorField{T,N}) where{T:<Number,N}

      Returns the tensor field dimension: N-1

"""
function tensortype(f::AbstractTensorField{T,N}) where{T<:Number,N}
  return T
end

#----------------------------------------------------------------------       


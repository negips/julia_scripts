#---------------------------------------------------------------------- 
@doc raw"""

        function tensorOP3D!(To::AbstractArray{T,N},Ti::AbstractArray{T,N},OPX::AbstractMatrix{T},OPY::AbstractMatrix{T},OPZ::AbstractMatrix{T},W::AbstractArray{T,N}) where {T<:Number,N}
          
        Operation (OPZ ⨂ OPY ⨂ OPX) on three-dimensional tensor arrays (Ti).

"""
@views function tensorOP3D!(To::AbstractArray{T,N},Ti::AbstractArray{T,N},OPX::AbstractMatrix{T},OPY::AbstractMatrix{T},OPZ::AbstractMatrix{T},W::AbstractArray{T,N}) where {T<:Number,N}

  si  = size(Ti)
  so  = size(To)
  sw  = size(W)
  sx  = size(OPX)
  sy  = size(OPY)
  sz  = size(OPZ)
 
  vTo = view(To,1:sx[1],1:si[2],1:si[3])
  tensor1OP!(vTo,Ti,OPX)

  vW  = view(W ,1:sx[1],1:sy[1],1:si[3])
  vTo = view(To,1:sx[1],1:si[2],1:si[3])
  tensor2OP!(vW ,vTo,OPY)

  vW  = view(W ,1:sx[1],1:sy[1],1:si[3])
  vTo = view(To,1:sx[1],1:sy[1],1:sz[1])
  tensor3OP!(vTo,vW ,OPZ)

  return nothing
end
#---------------------------------------------------------------------- 
@doc raw"""

        function tensorOP2D!(To::AbstractArray{T,N},Ti::AbstractArray{T,N},OPX::AbstractMatrix{T},OPY::AbstractMatrix{T},W::AbstractArray{T,N}) where {T<:Number,N}
          
        Operation (OPY ⨂ OPX) on two-dimensional tensor arrays (Ti).

"""
@views function tensorOP2D!(To::AbstractArray{T,N},Ti::AbstractArray{T,N},OPX::AbstractMatrix{T},OPY::AbstractMatrix{T},W::AbstractArray{T,N}) where {T<:Number,N}

  si  = size(Ti)
  so  = size(To)
  sw  = size(W)
  sx  = size(OPX)
  sy  = size(OPY)
 
  vW  = view(W,1:sx[1],1:si[2])
  tensor1OP!(vW,Ti,OPX)

  vTo  = view(To,1:sx[1],1:sy[1])
  vW   = view(W ,1:sx[1],1:si[2])
  tensor2OP!(vTo,vW,OPY)

  return nothing
end
#---------------------------------------------------------------------- 

@doc raw"""

        function tensor1OP!(To::AbstractArray{T,N},Ti::AbstractArray{T,N},OP::AbstractMatrix{T}) where {T<:Number,N}
          
        Operation (OP) on multidimensional tensor arrays (Ti).
        OP is applied along the 1st dimension.

"""
@views function tensor1OP!(To::AbstractArray{T,N},Ti::AbstractArray{T,N},OP::AbstractMatrix{T}) where {T<:Number,N}

  si  = size(Ti)
  so  = size(To)
  r,c = size(OP)

  @assert si[1] == c "Incompatible Tensor and Operator dimensions." 
  @assert so[1] == r "Incompatible Tensor and Operator dimensions." 

  if (N == 1)
    # To = OP*Ti
    mul!(To,OP,Ti)
  else
    p   = prod(si[2:N])
   
    tup = Tuple(si[2:N])
    CI  = CartesianIndices(tup)
    for ci in CI
      mul!(To[:,ci],OP,Ti[:,ci])
    end  
  end  

  return nothing
end
#----------------------------------------------------------------------       
@doc raw"""

        function tensor2OP!(To::AbstractArray{T,N},Ti::AbstractArray{T,N},OP::AbstractMatrix{T}) where {T<:Number,N}
          
        Operation (OP) on multidimensional tensor arrays (Ti).
        OP is applied along the 1st dimension.

"""
@views function tensor2OP!(To::AbstractArray{T,N},Ti::AbstractArray{T,N},OP::AbstractMatrix{T}) where {T<:Number,N}

  si  = size(Ti)
  so  = size(To)
  r,c = size(OP)

  # @assert si[1] == c "Incompatible Tensor and Operator dimensions." 
  # @assert so[1] == r "Incompatible Tensor and Operator dimensions." 

  # To = Ti*OP'
  if ( N > 2)
    tup = Tuple(si[3:N])
    CI  = CartesianIndices(tup)
    OPT = OP'
  
    for ci in CI
      mul!(To[:,:,ci],Ti[:,:,ci],OPT)
    end
  else
    OPT = OP'
    mul!(To,Ti,OPT)
  end    

  return nothing
end
#----------------------------------------------------------------------       
@doc raw"""

        function tensor3OP!(To::AbstractArray{T,N},Ti::AbstractArray{T,N},OP::AbstractMatrix{T}) where {T<:Number,N}
          
        Operation (OP) on multidimensional tensor arrays (Ti).
        OP is applied along the 1st dimension.

"""
@views function tensor3OP!(To::AbstractArray{T,N},Ti::AbstractArray{T,N},OP::AbstractMatrix{T}) where {T<:Number,N}

  si  = size(Ti)
  so  = size(To)
  r,c = size(OP)

  # @assert si[1] == c "Incompatible Tensor and Operator dimensions." 
  # @assert so[1] == r "Incompatible Tensor and Operator dimensions." 

  # To = Ti*OP'
  if ( N > 3)
    p12 = si[1]*si[2]
    RTi = reshape(Ti,p12,si[3:N])
 
    p12 = so[1]*so[2]
    RTo = reshape(To,p12,so[3:N])

    tup = Tuple(si[4:N])
    CI  = CartesianIndices(tup)
    OPT = OP'

    for ci in CI
      mul!(RTo[:,:,ci],RTi[:,:,ci],OPT)
    end
  else

    p12 = si[1]*si[2]
    RTi = reshape(Ti,p12,si[3])
 
    p12 = so[1]*so[2]
    RTo = reshape(To,p12,so[3])

    OPT = OP'
    mul!(RTo,RTi,OPT)

  end    

  return nothing
end
#----------------------------------------------------------------------       





# @doc raw"""
# 
#         function tensor2OP!(To::AbstractArray{T,N},Ti::AbstractArray{T,N},OP::AbstractMatrix{T},W::AbstractArray{T,N}) where {T<:Number,N}
#           
#         Operation (OP) on two-dimensional tensor arrays (tfin).
#         OP is applied along all dimensions.
# 
# """
# @views function tensor2OP!(To::AbstractArray{T,N},Ti::AbstractArray{T,N},OP::AbstractMatrix{T},W::AbstractArray{T,N}) where {T<:Number,N}
# 
#   ri,ci    = size(Ti)
#   ro,co    = size(To)
#   r,c      = size(OP)
# 
#   w1       = W[1:r,1:ci]
#   # W = OP*Ti
#   mul!(w1,OP,Ti)
# 
#   # To = W*(OP^T)
#   OPT = OP'
#   mul!(To,w1,OPT)
# 
#   return nothing
# end
# #----------------------------------------------------------------------       
# 
# @doc raw"""
# 
#         function tensor3OP!(tfout::AbstractArray{T,N},tfin::AbstractArray{T,N},OP::AbstractMatrix{T},wk::AbstractArray{T,N}) where {T<:Number}
#           
#         Operation (OP) on three-dimensional tensor arrays (tfin).
#         OP is applied along all dimensions.
# 
# """
# @views function tensor3OP!(To::AbstractArray{T,N},Ti::AbstractArray{T,N},OP::AbstractMatrix{T},W::AbstractArray{T,N}) where {T<:Number,N}
# 
#          si       = size(Ti)
#          so       = size(To)
#          sop      = size(OP)
# 
#          w1       = To[1:sop[1],1:si[2],1:si[3]]   # slices are views
#          tup      = (sop[1],si[2]*si[3])
#          w2       = reshape(w1,tup)
# 
#          tup      = (si[1],si[2]*si[3])
#          w3       = reshape(Ti,tup)
# 
#          # 1st Dimension
#          # W = OP*Ti
#          mul!(w2,OP,w3)
# 
#          # 2nd Dimension
#          OPT = OP'
#          for k in 1:si[3]
# 
#            v4 = wk[1:sop[1],1:sop[1],k]
#            v5 = v1[1:sop[1],1:si[2],k]
#            mul!(v4,v5,OPT) 
# 
#          end
#          
#          # 3rd Dimension
#          v6  = wk[1:sop[1],1:sop[1],1:si[3]]
#          tup = (sop[1]*sop[1],si[3])
#          v7  = reshape(v6,tup)
#          tup = (sop[1]*sop[1],sop[1])
#          v8  = reshape(tfout,tup)
#          mul!(v8,v7,OPT)
# 
#          return nothing
#        end
# #----------------------------------------------------------------------       


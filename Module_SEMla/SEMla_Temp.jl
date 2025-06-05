#     Non-Constructor Function Definitions
#---------------------------------------------------------------------- 
@doc raw"""

        function tensor2OP!(tfout::AbstractArray{T,N},tfin::AbstractArray{T,N},OP::AbstractMatrix{T},wk::AbstractArray{T,N}) where {T<:Number,N}
          
        Operation (OP) on two-dimensional tensor arrays (tfin).
        OP is applied along all dimensions.

"""
@views function tensor2OP!(tfout::AbstractArray{T,N},tfin::AbstractArray{T,N},OP::AbstractMatrix{T},wk::AbstractArray{T,N}) where {T<:Number,N}

  si       = size(tfin)
  so       = size(tfout)
  sop      = size(OP)

  v1       = wk[1:sop[1],1:si[2]]
  # wk = OP*tfin
  mul!(v1,OP,tfin)

  # tfout = wk*(OP^T)
  OPT = OP'
  mul!(tfout,v1,OPT)

  return nothing
end
#----------------------------------------------------------------------       

@doc raw"""

        function tensor3OP!(tfout::AbstractArray{T,N},tfin::AbstractArray{T,N},OP::AbstractMatrix{T},wk::AbstractArray{T,N}) where {T<:Number}
          
        Operation (OP) on three-dimensional tensor arrays (tfin).
        OP is applied along all dimensions.

"""
@views function tensor3OP!(tfout::AbstractArray{T,N},tfin::AbstractArray{T,N},OP::AbstractMatrix{T},wk::AbstractArray{T,N}) where {T<:Number,N}

         si       = size(tfin)
         so       = size(tfout)
         sop      = size(OP)

         v1       = tfout[1:sop[1],1:si[2],1:si[3]]   # slices are views
         tup      = (sop[1],si[2]*si[3])
         v2       = reshape(v1,tup)

         tup      = (si[1],si[2]*si[3])
         v3       = reshape(tfin,tup)

         # 1st Dimension
         mul!(v2,OP,v3)

         # 2nd Dimension
         OPT = OP'
         for k in 1:si[3]

           v4 = wk[1:sop[1],1:sop[1],k]
           v5 = v1[1:sop[1],1:si[2],k]
           mul!(v4,v5,OPT) 

         end
         
         # 3rd Dimension
         v6  = wk[1:sop[1],1:sop[1],1:si[3]]
         tup = (sop[1]*sop[1],si[3])
         v7  = reshape(v6,tup)
         tup = (sop[1]*sop[1],sop[1])
         v8  = reshape(tfout,tup)
         mul!(v8,v7,OPT)

         return nothing
       end
#----------------------------------------------------------------------       
@views function tensor3OP_2!(tfout::AbstractArray{T,N},tfin::AbstractArray{T,N},OP::AbstractMatrix{T},wk::AbstractArray{T,N}) where {T<:Number,N}

           si       = size(tfin)
           so       = size(tfout)
           sop      = size(OP)

           u0       = view(tfout,:,:,:)
           v0       = view(tfin,:,:,:)
           w0       = view(wk,:,:,:)

           # Slices are views 
           u1       = u0[1:sop[1],1:si[2],1:si[3]]
           tup      = (sop[1],si[2]*si[3])
           u2       = reshape(u1,tup)

           tup      = (si[1],si[2]*si[3])
           v1       = reshape(v0,tup)

           # 1st Dimension
           mul!(u2,OP,v1)

           # 2nd Dimension
           OPT = OP'
           for k in 1:si[3]

             w1 = w0[1:sop[1],1:sop[1],k]
             u3 = u0[1:sop[1],1:si[2],k]
             mul!(w1,u3,OPT) 

           end
           
           # 3rd Dimension
           w2  = w0[1:sop[1],1:sop[1],1:si[3]]
           tup = (sop[1]*sop[1],si[3])
           w3  = reshape(w2,tup)

           tup = (sop[1]*sop[1],sop[1])
           u4  = reshape(u0,tup)
           mul!(u4,w3,OPT)

         return nothing
       end
#----------------------------------------------------------------------       

#----------------------------------------------------------------------       






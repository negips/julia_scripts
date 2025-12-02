@doc """
      NLGinzburgLandau(L::AbstractMatrix{T},B::AbstractMatrix{T},v::AbstractVector{T},δ::T,bc1::T,bc2::T,ifbc1::Bool,ifbc2::Bool)

      Calculate Fully Non-linear time derivative operator (with mass matrix).

"""
function NLGinzburgLandau(L,B,v,δ::T,bc1::T,bc2::T,ifbc1::Bool,ifbc2::Bool) where {T <: Number}

   dv  = L*v
   for i in LinearIndices(dv)
     dv[i] = dv[i] + δ*B[i]*conj(v[i])*v[i]*v[i]
   end

   # Note these are boundary conditions for time derivative.
   # Not for the field values themselves.
   if (ifbc1) 
     dv[1] = bc1
   end

   if (ifbc2) 
     dv[end] = bc2
   end  

   return dv
end  
#---------------------------------------------------------------------- 
@doc """
      ForcedNLGinzburgLandau(L,B,v::AbstractVector{T},F,θ,δ::T,ω,bc1::T,bc2::T,ifbc1::Bool,ifbc2::Bool)

      Calculate Fully Non-linear harmonically forced time derivative operator (with mass matrix).

"""
function ForcedNLGinzburgLandau(L::AbstractMatrix{T},B,v::AbstractVector{T},F::AbstractMatrix{T},θ::AbstractVector{T},δ,Ω::AbstractVector{T},bc1,bc2,ifbc1::Bool,ifbc2::Bool) where {T <: Number}

   dv = L*v
   fθ = F*θ
   for i in LinearIndices(dv)
     dv[i] = dv[i] + δ*B[i]*conj(v[i])*v[i]*v[i] + fθ[i]
   end
   # Note these are boundary conditions for time derivative.
   # Not for the field values themselves.
   if (ifbc1) 
     dv[1] = bc1
   end
   if (ifbc2) 
     dv[end] = bc2
   end  

   dθ = Ω.*θ 

   return dv, dθ
end  
#---------------------------------------------------------------------- 








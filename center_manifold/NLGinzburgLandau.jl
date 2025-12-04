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
@doc """
StuartLandau3(G1::AbstractMatrix{T},G2::AbstractMatrix{T},G3::AbstractMatrix{T},v::AbstractVector{T}) where {T <: Number}

      Calculate the time derivative of the (third-order) Stuart-Landau equations.

"""
function StuartLandau3(G1::AbstractMatrix{T},G2::AbstractMatrix{T},G3::AbstractMatrix{T},v::AbstractVector{T}) where {T <: Number}

   nv,nt1 = size(G1)
   ord1   = 1
   dv     = G1*v

   nv,nt2 = size(G2)
   ord1   = 1

   v2     = zeros(T,nt2)
   for i in 1:nt2



   end        


   return dv
end  
#---------------------------------------------------------------------- 








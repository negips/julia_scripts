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
     dv[i] = dv[i] + δ*B[i]*conj(v[i])*v[i]*v[i] + B[i]*fθ[i]
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
      DampedFGL(OPg,B,v::AbstractVector{T},V::AbstractMatrix{T},W::AbstractMatrix{T},χ::AbstractVector{T})

      Calculate Damped, Fully Non-linear, harmonically forced time derivative operator (with mass matrix).

"""
function DampedFGL(FGL,B::AbstractVector{S},v::AbstractVector{T},θ::AbstractVector{T},V::AbstractMatrix{T},W::AbstractMatrix{T},χ::AbstractVector{T}) where {T,S <: Number}

   dv,dθ = FGL(v,θ)

   # Projections
   α     = W'*(B.*v)
   dv   .= dv .- B.*(V*(α.*χ))

   return dv, dθ
end  
#---------------------------------------------------------------------- 

@doc """
      StuartLandau1(G1::AbstractMatrix{T},v::AbstractVector{T}) where {T <: Number}

      Calculate the time derivative of the (first-order) Stuart-Landau equations.

"""
function StuartLandau1(G1::AbstractMatrix{T},v::AbstractVector{T}) where {T <: Number}

   ord1     = 1
   dv       = G1*v

   return dv
end  
#---------------------------------------------------------------------- 
@doc """
      StuartLandau3(G1::AbstractMatrix{T},G2::AbstractMatrix{T},G3::AbstractMatrix{T},v::AbstractVector{T}) where {T <: Number}

      Calculate the time derivative of the (third-order) Stuart-Landau equations.

"""
function StuartLandau3(G1::AbstractMatrix{T},G2::AbstractMatrix{T},G3::AbstractMatrix{T},v::AbstractVector{T}) where {T <: Number}

   ord1     = 1
   dv       = G1*v

   Ord      = 2
   NLv2     = CenterManifold.EvaluateNonLinear(v,Ord)
   dv      .= dv .+ G2*NLv2

   Ord      = 3
   NLv3     = CenterManifold.EvaluateNonLinear(v,Ord)
   dv      .= dv .+ G3*NLv3

   return dv
end  
#---------------------------------------------------------------------- 
@doc """
      StuartLandau5(G1::AbstractMatrix{T},G2::AbstractMatrix{T},G3::AbstractMatrix{T},G4::AbstractMatrix{T},G5::AbstractMatrix{T},v::AbstractVector{T}) where {T <: Number}

      Calculate the time derivative of the (third-order) Stuart-Landau equations.

"""
function StuartLandau5(G1::AbstractMatrix{T},G2::AbstractMatrix{T},G3::AbstractMatrix{T},G4::AbstractMatrix{T},G5::AbstractMatrix{T},v::AbstractVector{T}) where {T <: Number}

   ord1     = 1
   dv       = G1*v

   Ord      = 2
   NLv2     = CenterManifold.EvaluateNonLinear(v,Ord)
   dv      .= dv .+ G2*NLv2

   Ord      = 3
   NLv3     = CenterManifold.EvaluateNonLinear(v,Ord)
   dv      .= dv .+ G3*NLv3

   Ord      = 4
   NLv4     = CenterManifold.EvaluateNonLinear(v,Ord)
   dv      .= dv .+ G4*NLv4

   Ord      = 5
   NLv5     = CenterManifold.EvaluateNonLinear(v,Ord)
   dv      .= dv .+ G5*NLv5

   return dv
end  
#---------------------------------------------------------------------- 








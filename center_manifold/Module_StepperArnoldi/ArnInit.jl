function ArnKrylovInit(StpInp::StepperInput,ArnInp::ArnoldiInput;Dtype=Float64)

      # Prec        = typeof(ArnInp).parameters[1]
      ndof        = ArnInp.vlen
      Nev         = ArnInp.nev                  # Number of eigenvalues to calculate
      EKryl       = ArnInp.ekryl                # Additional size of Krylov space
      LKryl       = ArnInp.lkryl                # Total Size of Krylov space    
      
      V           = zeros(Dtype,ndof,LKryl+1)
      H           = zeros(Dtype,LKryl+1,LKryl)
      
      return V,H
end      
#---------------------------------------------------------------------- 
function StpArn_SetBC!(v::AbstractVector{T},lbc::Bool,rbc::Bool) where {T<:Number}

  if lbc
    v[1]   = T(0)
  end
  if rbc
    v[end] = T(0)
  end

  return nothing
end
#---------------------------------------------------------------------- 
function ArnInitVector(vlen::Int,lbc::Bool,rbc::Bool;Dtype=Float64)

  rng       = Xoshiro(123)
  v         = rand(rng,Dtype,vlen)

  StpArn_SetBC!(v,lbc,rbc)

  return v
end  





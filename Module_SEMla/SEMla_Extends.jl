# This is where we add extensions of Julia functions
#---------------------------------------------------------------------- 
function copy(f::TensorField{T,N}) where {T<:Number,N}

  tfield    = Base.copy(f.tfield)

  return TensorField(tfield)
end
#----------------------------------------------------------------------
function copy(vf::Vector{TensorField{T,N}}) where {T<:Number,N}

  nt        = length(vf)
  tfields   = Vector{TensorField{T,N}}(undef,nt)
  for i in 1:nt
    tfields[i] = copy(vf[i])
  end  

  return tfields
end        
#----------------------------------------------------------------------       
function copy(f::NTensorFields{T,N}) where {T<:Number,N}

  n         = f.ntensors
  TF        = copy(f.TFields)

  return NTensorFields(TF)
end
#----------------------------------------------------------------------       





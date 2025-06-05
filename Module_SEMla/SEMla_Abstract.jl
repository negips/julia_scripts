# Add Definition of new structures here
abstract type AbstractTensorField{T,N} <: AbstractArray{T,N} end
abstract type AbstractIsoTensorField{T,N} <: AbstractTensorField{T,N} end
#---------------------------------------------------------------------- 

# Specialized Type aliases for convenience
"""
    AbstractTensorField2{T}

Supertype for two-dimensional AbstractTensorFields with
elements of type `T`. Alias for [`AbstractTensorField{T,3}`](@ref).
"""
const AbstractTensorField2{T} = AbstractTensorField{T,3}
#---------------------------------------------------------------------- 
"""
    AbstractTensorField3{T}

Supertype for two-dimensional AbstractTensorFields with
elements of type `T`. Alias for [`AbstractTensorField{T,3}`](@ref).
"""
const AbstractTensorField3{T} = AbstractTensorField{T,4}
#---------------------------------------------------------------------- 

#      export AbstractTensorfield,
#             AbstractIsoParTensorField

#      abstract type TwoTensorField <: NTensorField end
#----------------------------------------------------------------------  

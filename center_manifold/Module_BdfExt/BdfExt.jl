# I should turn this into a module when i'm done editing

module BdfExt

  # using LinearAlgebra
  # using Printf

  include("GetBDF.jl")
  include("BdfExtStruct.jl")
  # include("BdfExtFunctions.jl")

  # Defined Structures
  export BdFields,
         ExtFields

  # Defined unctions
  export GetBDF,
         GetBDFT2,
         GetEXT

  # Dynamical Evolution related functions 
  # export DynamicalSystem1


end
#---------------------------------------------------------------------- 











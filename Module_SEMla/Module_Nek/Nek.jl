#     Author:     Prabal Negi
#---------------------------------------------------------------------- 
      module Nek

      import Base.copy
      import Base.copy!
      import Base.copyto!

      using LinearAlgebra
      using PolynomialBases
      using MPI
      using HDF5
      using Printf
      using Statistics        # Refinements

#     Basic definitions: 
#     abstract types, structures, constructors, extensions
#--------------------------------------------------      
      include("Nek_Structs.jl")
      include("Nek_Constructors.jl")
      include("Nek_Extends.jl")
      include("Nek_ElementRefine.jl")

#     Define Functions
#--------------------------------------------------      
      include("Nek_ParamDescriptions.jl")
      include("Nek_Functions.jl")
      include("Nek_WriteRea.jl")

      # Nek Based Structures/Functions
#--------------------------------------------------     

      # Structures/Constructors 
      export NekField,
             Re2Field,
             Ma2Field,
             NekObjects,
             ReaData

      # Functions      
      export read_re2,                    # Multiple Definitions
             read_ma2,                    # Multiple Definitions
             read_fld,                    # Multiple Definitions
             param_descriptions,          # Discriptions of Parameters
             write_reafile,
             element_refine2D,            # Refine 2D element
             twoelements_refine2D         # Refine 2 elements (2D).


#----------------------------------------------------------------------  
      function __init__()
      # Initialization at Module load        

        return nothing
      end 

#----------------------------------------------------------------------


      function nek_init()
      # Manual initialization 

        return nothing
      end 

#----------------------------------------------------------------------


#---------------------------------------------------------------------- 
      end   # Module Nek











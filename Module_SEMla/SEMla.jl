#     Author:     Prabal Negi
#

      module SEMla

      import Base.copy
      import Base.copy!
      import Base.copyto!

      using LinearAlgebra
      using PolynomialBases
      using MPI
      using HDF5
#      using UnsafeArrays

      include("Module_Nek/Nek.jl")
      using .Nek

#     Basic definitions: 
#     abstract types, structures, constructors, extensions
#--------------------------------------------------      
      include("SEMla_Abstract.jl")        # 
      include("SEMla_Structs.jl")         # 
      include("SEMla_Constructors.jl")    #
      include("SEMla_Extends.jl")         # native function extensions
      include("SEMla_TensorFunctions.jl") # Tensor-functions
      include("SEMla_Functions.jl")

      include("SEMla_Temp.jl")            # Temporary functions before they are moved to the main module

##     Maybe I should make this another module?
#      include("Nek_Structs.jl")
#      include("Nek_Constructors.jl")
#      include("Nek_Extends.jl")

#     Define Functions
#--------------------------------------------------      
#      include("SEMla_Nek.jl")
      include("SEMla_Mesh.jl")

      # Abstraact Definitions      
      export AbstractTensorField,
             AbstractIsoParTensorField,
             AbstractTensorField2,
             AbstractTensorField3
            

      # Structures (and Constructors) 
      export SEMla_Config,
             TensorField,
             IsoTensorField,
             NTensorFields,
             LocalVertexMap,
             GlobalVertexMap1S,
             GlobalVertexMap,
             GlobalCommMap,
             LocalEdgeMap,
             GlobalEdgeMap

      # Functions
      export semla_init,
             gen_rema2,
             distribute_mesh,
             distribute_unsorted_mesh,
             tensorOPn!,
             tensorfieldOP!,
             tensorfieldOP12!,
             tensordimension,
             localuniquemap,
             localuniquetoglobal,
             buildlocalvertexmap,
             buildglobalvertexmap,
             buildlocaledgemap,
             buildglobaledgemap,
             edgeindices



      # Nek Based Structures/Functions
#--------------------------------------------------     

      # Structures      
      export NekField,
             Re2Field,
             Ma2Field

      # Functions      
      export read_re2,        # Multiple Definitions
             read_ma2,        # Multiple Definitions
             read_fld         # Multiple Definitions

#----------------------------------------------------------------------  
      function __init__()
      # Initialization at Module load        

        return nothing
      end 

#----------------------------------------------------------------------


      function semla_init()
      # Manual initialization 

        return nothing
      end 

#----------------------------------------------------------------------


#---------------------------------------------------------------------- 
      end   # Module SEMla











#     Author:     Prabal Negi
#

      module SEM1D

      import Base.copy
      import Base.copy!
      import Base.copyto!

      using LinearAlgebra
      using PolynomialBases
      using SparseArrays


#     Basic definitions: 
#     abstract types, structures, constructors, extensions
#--------------------------------------------------      
      include("SEM1D_Abstract.jl")        # 
      include("SEM1D_Structs.jl")         # 
      include("SEM1D_Constructors.jl")    #
      include("SEM1D_Geom.jl")            #


#    Abstracts  
#-------------------------------------------------- 
      export AbstractParams,
             AbstractGeometry

#     Structures (and Constructors) 
#--------------------------------------------------      

      export SEM_Input,
             SEM_GeoMat


#     Functions
#--------------------------------------------------      


      end         # Module SEM1D       
#---------------------------------------------------------------------- 





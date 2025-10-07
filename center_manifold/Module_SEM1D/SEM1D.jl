#     Author:     Prabal Negi
#

      module SEM1D

      import Base.copy
      import Base.copy!
      import Base.copyto!

      using LinearAlgebra
      using PolynomialBases
      using SparseArrays
      using SpecialFunctions
      using Roots


#     Basic definitions: 
#     abstract types, structures, constructors, extensions
#--------------------------------------------------      
      include("SEM1D_Abstract.jl")        # 
      include("SEM1D_Structs.jl")         # 
      include("SEM1D_Constructors.jl")    #
      include("SEM1D_Geom.jl")            #
      include("GinzburgLandau.jl")        #
      include("SEM1D_QQT.jl")             #


#    Abstracts  
#-------------------------------------------------- 
      export AbstractParams,
             AbstractGeometry

#     Structures (and Constructors) 
#--------------------------------------------------      

      export SEMInput,
             SEMGeoMat


#     Functions
#--------------------------------------------------      

      export GinzburgLandauSparse,
             AdjointGinzburgLandauSparse,
             GLSetBC!,
             SEM_SetBC!,
             GLAnalyticalSpectra,
             SEM_QQT,
             SEM_Global_Num

      end         # Module SEM1D       
#---------------------------------------------------------------------- 





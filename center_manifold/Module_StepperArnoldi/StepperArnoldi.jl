      module StepperArnoldi

      import Base.copy
      # import Base.copy!
      # import Base.copyto!

      using LinearAlgebra
      using IterativeSolvers
      using Roots
      using Random
      using JLD2
      using Printf
      
#     Basic definitions: 
#     abstract types, structures, constructors, extensions
#--------------------------------------------------
      include("StpArn_Abstract.jl")        # 
      include("StpArn_Struct.jl")          # 

      include("StpArnFunctions.jl")
      include("RK4.jl")
      include("ArnUpd.jl")
      include("ArnIRst.jl")
      include("ExplicitShiftedQR.jl")
      include("BulgeChase.jl")
      include("IRAM.jl")
      include("ArnInit.jl")
      include("Stepper.jl")


#      include("SEM1D_Constructors.jl")    #
#      include("SEM1D_Geom.jl")            #
#      include("GinzburgLandau.jl")        #
#      include("SEM1D_QQT.jl")             #

#     Abstracts  
#-------------------------------------------------- 
      export AbstractStepperParams,
             AbstractArnoldiParams

#     Structures (and Constructors) 
#--------------------------------------------------      
      export StepperInput,
             ArnoldiInput,
             ArnoldiOutput


#     Functions
#--------------------------------------------------      
      export RK4!,
             BiRK4!,
             RestrictedBRK4!,
             ArnUpd,
             ArnIRst,
             ArnBIRst,
             ArnGetUpperShifts,
             ArnGetLowerShifts,
             ArnGetCustomShifts,
             ExplicitShiftedQR,
             FrancisAlg,
             FrancisSeq,
             FrancisSeqExact,
             RevFrancisSeq,
             CreateBulge,
             CreateLowerBulge,
             polyHe1,
             enpolyH,
             polyAe1,
             CreateReflector,
             CreateReflectorZeros,
             CreateGivens,
             CreateReflectorZeros2,
             TransposeReflectorZeros,
             AdjointReflectorZeros,
             BandedHessenberg,
             ChaseBulgeDown,
             ChaseBulgeDown1,
             ChaseBulgeUp,
             IRAM!

      export ArnKrylovInit,
             StpArn_SetBC!,
             ArnInitVector,
             StepArn,
             ObliqueSubspaceRemoval!,
             RestrictedStepArn


      end         # Module StepperArnoldi
#----------------------------------------------------------------------




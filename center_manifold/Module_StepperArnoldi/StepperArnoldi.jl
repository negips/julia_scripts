      module StepperArnoldi

      import Base.copy
      # import Base.copy!
      # import Base.copyto!

      using SparseArrays
      using LinearAlgebra
      using IterativeSolvers
      using Roots
      using Random
      using JLD2
      using Printf
      
#     Basic definitions: 
#     abstract types, structures, constructors, extensions
#--------------------------------------------------
      include("StpArnAbstract.jl")        # 
      include("StpArnStruct.jl")          # 
      include("StpArnConstructors.jl")    # 

      include("StpArnFunctions.jl")
      include("RK4.jl")
      include("OPRK4.jl")
      include("ArnUpd.jl")
      include("ArnIRst.jl")
      include("ExplicitShiftedQR.jl")
      include("BulgeChase.jl")
      include("IRAM.jl")
      include("ArnInit.jl")
      include("Stepper.jl")
      include("StepperOP.jl")
      include("ExtendTangentSpace.jl")
      include("ExtendPertTangentSpace.jl")


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
             ArnoldiOutput,
             ExtendedMode,
             ExtendedModes


#     Functions
#--------------------------------------------------      
      export RK4!,
             BiRK4!,
             RestrictedBRK4!,
             OPRK4!,
             OPBiRK4!,
             OPRestrictedBRK4!,
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
             StepArnOP,
             RestrictedStepArn,
             RestrictedStepArnOP,
             REPStepArn,
             ObliqueSubspaceRemoval!,
             ObliqueSubspaceRemoval2!,
             ObliqueSubspaceRemoval3!,
             PLx,
             PLx!,
             RPLx!,
             EPLx,
             REPLx,
             ExtLOP,
             PertLx,
             PertLx!,
             ExtPertLx,
             RestrictedExtPertLx,
             ExtendTangentSpace,
             ExtendTangentSpaceOP,
             ExtendTangentSpaceRestricted,
             ExtendTangentSpaceRestrictedOP,
             ExtendedTangentSpaces,
             ExtendedTangentSpacesOP,
             EPTangentSpaces,
             EPTangentSpace

      end         # Module StepperArnoldi
#----------------------------------------------------------------------




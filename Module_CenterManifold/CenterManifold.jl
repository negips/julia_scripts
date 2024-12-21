# I should turn this into a module when i'm done editing

module CenterManifold

  using SymPy
  using LinearAlgebra

  # These two are for reading NekFiles
  # We could possibly move the Nek related functions out.
  using CSV         
  using DataFrames

  include("CMStruct.jl")
  include("CMFunctions.jl")
  include("CMSymbolic.jl")
  include("MapSymbols.jl")
  include("CMDynamicalSystem.jl")

  # Defined Structures
  export ROPMap

  # Non-symbolic functions
  export GetNCritical,
         ReadNekCMFiles,
         BuildReducedOPMap,
         FillPowers!,
         GetPowers,
         GetNonZeroPowers,
         QuadTuples,
         InteractionNTuples,
         NInteractionTerms,
         NInteractionTermsNZ,
         TotalInteractionTerms,
         NextRow,
         NextRow!,
         UpdatePolynomialIndex!,
         newindex,
         oldindex,
         GetInteractionTerms,
         GetInteractionIndex

  # Dynamical Evolution related functions 
  export DynamicalSystem1,
         EvolveSystem1,
         DynamicalSystemN,
         EvolveSystemN!
        

  # Symbolic functions
  export CollectTerms!,
         CollectTerms,
         CollectOrderN!,
         GenerateMR,
         GetSolnExpressions,
         GetCoeffs,
         GetPQMap,
         BuildGenericEquation,
         AddLinearOP!,
         GetNonLinearTerms,
         GetAllNonLinearTerms,
         GenerateSymbols,
         GenerateSymbolsSuper,
         GenerateSymMatrix,
         GenerateNonLinearSymbols,
         GenerateAllNonLinearSymbols,
         GenDict0,
         GenDict,
         SymbolFunction

  # Mapping of symbols
  export Map2Values,
         MapQ2ABC,
         Mapθ2λ,
         MapX2CSP,
         MapSolutions


end
#---------------------------------------------------------------------- 











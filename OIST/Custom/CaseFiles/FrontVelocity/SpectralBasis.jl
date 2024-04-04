module SpectralBases

  using PolynomialBases: AbstractDomain, Line, AbstractBasis, laguerre
  using FastGaussQuadrature: gausslaguerre
  using AutoHashEquals: @auto_hash_equals
  using ArgCheck

  abstract type SpectralBasis{Domain} <: AbstractBasis{Domain} end

  include("Laguerre.jl")

  export SpectralBasis

  export LaguerreSpectral,LaguerreDerivative,LaguerreDerivativeMatrix,
         LaguerreTransforms

end         # module









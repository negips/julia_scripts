#     Testing Chebyshev differentiation
#

      using LinearAlgebra,ToeplitzMatrices,PyPlot

      include("chebdiff.jl")
      include("SoarMatrix.jl")
      include("SoarMatrixOS.jl")

      close()

#      N=5;
#      M=4;
#
#      x,DM = chebdif(N,M);
#
#      plot(x)
#

      N = 5;
      ω = 0.3;
      β = 1.0;
#      x,y,U,DM,K,D,M,k1,k2,k3 = SoarMatrix(ω,N)
      x,y,U,DM,K,D,M,k1,k2 = SoarMatrixOS(ω,β,N)

      plot(y,U,marker=".")


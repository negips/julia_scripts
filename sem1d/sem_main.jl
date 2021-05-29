println("Main interface for 1D SEM")

using PolynomialBases
using LinearAlgebra
# using UnicodePlots
using Plots

# Include the function files
include("sem_geom.jl")

include("AssembleMatrix.jl")

include("AssembleMatrixLesshafft.jl")


include("sem_init_ref.jl")

c0 = 1.0;
ν  = 0.1;   # nu


#L  = AssembleMatrix(c0,Geom.cnv,ν,Geom.wlp,lx1,nel);

#L  = AssembleMatrixLesshafft(c0,Geom.cnv,Geom.wlp,Geom.xm1,Geom.bm1,Basis,lx1,nel);



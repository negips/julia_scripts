println("Main interface for 1D SEM")

using PolynomialBases

# Include the function files
include("sem_geom.jl")

include("sem_init_ref.jl")

include("AssembleMatrix.jl")

c0 = 1.0;
ν  = 0.1;

L  = AssembleMatrix(c0,Geom.cnv,ν,Geom.wlp,lx1,nel);





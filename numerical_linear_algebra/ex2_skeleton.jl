using Pkg
Pkg.add("MatrixDepot")
using MatrixDepot, Random
nn=10;
Random.seed!(0)
A=matrixdepot("wathen",nn,nn)

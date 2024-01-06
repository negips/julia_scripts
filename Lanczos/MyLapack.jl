module MyLapack

  using LinearAlgebra.BLAS: @blasfunc
  using LinearAlgebra: BlasInt
  using LinearAlgebra
  
  const libsrc = "/usr/lib/julia/libblastrampoline.so.5"    # doesn't work for some reason
  const localsrc = "liblapack"

  include("lapack_dlaexc.jl")

  export dlaexc!




end

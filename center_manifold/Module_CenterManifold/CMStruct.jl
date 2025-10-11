@doc raw"""
      mutable struct ROPMap

      Has fields:
      LOP::Matrix{ComplexF64}         # Real/Complex
      LOPr::Matrix{ComplexF64}        # Real/Complex
      λ::Vector{ComplexF64}           # Real/Complex
      V::Matrix{ComplexF64}           # Right Eigenvector Matrix
      Vs::Matrix                      # Symbolic Right Eigenvector Matrix
      W::Matrix{ComplexF64}           # Left Eigenvector Matrix
      Ws::Matrix                      # Symbolic Left Eigenvector Matrix
      MapV::Dict{Sym,ComplexF64}      # Using Vs/Ws symbols
      MapW::Dict{Sym,ComplexF64}      # Using Vs/Ws symbols
      X2YSym::Dict{Sym,Sym}           # Using Vs/Ws symbols
      X2Y::Dict{Sym,Sym}              # Using V/W values
      X2YSym::Dict{Sym,Sym}           # Using Vs/Ws symbols
      Y2X::Dict{Sym,Sym}              # Using V/W values
      Y2XSym::Dict{Sym,Sym}           # Using Vs/Ws symbols

"""
mutable struct ROPMap

  # SymPy.SymbolicObject <: Number

  LOP::Matrix{ComplexF64}           # Real/Complex
  LOPr::Matrix{ComplexF64}          # Real/Complex
  λ::Vector{ComplexF64}             # Real/Complex
  V::Matrix{ComplexF64}             # Right Eigenvector Matrix
  Vs::Matrix                        # Symbolic Right Eigenvector Matrix
  W::Matrix{ComplexF64}             # Left Eigenvector Matrix
  Ws::Matrix                        # Symbolic Left Eigenvector Matrix
  MapV::Dict{Sym,ComplexF64}        # Using Vs/Ws symbols
  MapW::Dict{Sym,ComplexF64}        # Using Vs/Ws symbols
  X2Y::Dict{Sym,Sym}                # Using V/W values
  X2YSym::Dict{Sym,Sym}             # Using Vs/Ws symbols
  Y2X::Dict{Sym,Sym}                # Using V/W values
  Y2XSym::Dict{Sym,Sym}             # Using Vs/Ws symbols

end
#---------------------------------------------------------------------- 


#     Add Definition of new structures here
#---------------------------------------------------------------------- 
"""
      mutable struct StepperInput

      Has fields:

      ifadjoint               - If Adjoint solve
      ifoptimal               - If Optimal perturbations
      ifverbose               - Verbose stepper
      verbosestep::Int        - Verbose step
      nsteps::Int             - Steps in one Stepper phase
      timestep::Float64       - DT


"""
mutable struct StepperInput <: AbstractStepperParams

      ifadjoint::Bool
      ifoptimal::Bool
      ifverbose::Bool
      verbosestep::Int 
      nsteps::Int
      timestep::Float64

end
#----------------------------------------------------------------------       
"""
      mutable struct ArnoldiInput

      Has fields:

      ifarnoldi         - If Enabled 
      ifverbose         - If Verbose
      vlen              - Vector Length
      nev               - Number of eigenvalues to calculate
      ekryl             - Additional size of Krylov space
      lkryl             - Total Size of Krylov space    
      eigshift          - Shift origin of eigenvalues
      ngs               - Number of Gram-Schmidt
      bsize             - Block Size
      outer_iterations  - Maximum no of outer loop iterations
      tol               - Convergence Tolerance 


"""
mutable struct ArnoldiInput <: AbstractStepperParams

      ifarnoldi::Bool 
      ifverbose::Bool
      ifeigshift::Bool
      vlen::Int
      nev::Int         
      ekryl::Int       
      lkryl::Int 
      eigshift::S where S<:Number
      ngs::Int
      bsize::Int
      outer_iterations::Int
      tol::Float64

end
#----------------------------------------------------------------------       
"""
      struct ArnoldiOutput

      Has fields:

      evals             - Eigenvalues
      evecs             - Eigenvectors
      ritz              - ritz values
      ifconv            - Convergence flag
      tol               - Convergence Tolerance 

"""
struct ArnoldiOutput{T}

      evals::Vector{T} 
      evecs::Matrix{T}
      ritz::Vector{T}
      ifconv::Bool  
      tol::Float64

end
#----------------------------------------------------------------------       
"""
      struct ExtendedMode

      Has fields:

      λe::T                   - Frequency of extended Mode.
      ve::AbstractVector{T}   - Direct Extended Mode.
      we::AbstractVector{T}   - Adjoint Extended Mode.
      z::AbstractMatrix{T}    - Original Adjoint Mode' extensions.
      Γ::AbstractVector{T}    - Jordan Form Couplings for resonant extensions.

"""
mutable struct ExtendedMode{T}

      λe::T
      fe::AbstractVector{T}
      Γ::AbstractVector{T}
      ve::AbstractVector{T}
      we::AbstractVector{T}
      z::AbstractMatrix{T}

end
#----------------------------------------------------------------------       
"""
      struct ExtendedModes

      (Plural)

      Has fields:

      λe::AbstractVector{T}   - Frequencies of Extended Modes.
      Γ::AbstractMatrix{T}    - Jordan Form Couplings for resonant extensions.
      Ve::AbstractMatrix{T}   - Direct Extended Modes.
      We::AbstractMatrix{T}   - Adjoint Extended Modes.
      Z::AbstractMatrix{T}    - Original Adjoint Modes' extensions.

"""
mutable struct ExtendedModes{T}

      λe::AbstractVector{T}
      Fe::AbstractMatrix{T}
      Γ::AbstractMatrix{T}
      Ve::AbstractMatrix{T}
      We::AbstractMatrix{T}
      Z::AbstractMatrix{T}

end
#----------------------------------------------------------------------       





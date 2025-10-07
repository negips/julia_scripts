#     Add Definition of new structures here
#---------------------------------------------------------------------- 
"""
      struct StepperInput

      Has fields:

      ifadjoint               - If Adjoint solve
      ifoptimal               - If Optimal perturbations
      ifverbose               - Verbose stepper
      verbosestep::Int        - Verbose step
      nsteps::Int             - Steps in one Stepper phase
      timestep::Float64       - DT


"""
struct StepperInput <: AbstractStepperParams

      ifadjoint::Bool
      ifoptimal::Bool
      ifverbose::Bool
      verbosestep::Int 
      nsteps::Int
      timestep::Float64

end
#----------------------------------------------------------------------       
"""
      struct ArnoldiInput

      Has fields:

      ifarnoldi         - If Enabled 
      ifverbose         - If Verbose
      vlen              - Vector Length
      nev               - Number of eigenvalues to calculate
      ekryl             - Additional size of Krylov space
      lkryl             - Total Size of Krylov space    
      ngs               - Number of Gram-Schmidt
      bsize             - Block Size
      outer_iterations  - Maximum no of outer loop iterations
      tol               - Convergence Tolerance 


"""
struct ArnoldiInput{T} <: AbstractStepperParams

      ifarnoldi::Bool 
      ifverbose::Bool
      vlen::Int
      nev::Int         
      ekryl::Int       
      lkryl::Int       
      ngs::Int
      bsize::Int
      outer_iterations::Int
      tol::T

end
#----------------------------------------------------------------------       
"""
      struct ArnoldiOutput

      Has fields:

      evals             - Eigenvalues
      evecs             - Eigenvectors
      ifconv            - Convergence flag
      tol               - Convergence Tolerance 


"""
struct ArnoldiOutput{T}

      evals::Vector{T} 
      evecs::Matrix{T}
      ifconv::Bool  
      tol::Float64

end
#----------------------------------------------------------------------       





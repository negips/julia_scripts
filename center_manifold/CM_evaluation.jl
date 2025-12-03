# Center Manifold Evaluation
println("Center-Manifold Evaluation on the Ginzburg Landau system.")


include("GL_Setup.jl")

include("GL_Standard_Tangent_Space.jl")

include("GL_Extended_Tangent_Space.jl")

include("GL_Asymptotic_CM_NormalForm.jl")

include("CM_WriteEquation.jl")

CM_DisplayTerms(1,"z",Khat,G_O2,G_O3)

println("Done.")








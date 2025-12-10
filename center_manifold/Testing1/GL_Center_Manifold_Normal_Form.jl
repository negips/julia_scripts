# Center Manifold Evaluation
println("Center-Manifold Evaluation on the Ginzburg Landau system.")

include("GL_Asymptotic_CM_NormalForm.jl")

CenterManifold.DisplayTerms3(1,"z",Khat,G_O2,G_O3)

println("Center-Manifold Evaluation Complete.")








# Center Manifold Evaluation
println("Center-Manifold Evaluation on the Ginzburg Landau system.")

ifnormal = false

if ifnormal
  include("GL_Asymptotic_CM_NormalForm.jl")
else
  include("GL_Asymptotic_CM_GraphForm.jl")
end

CenterManifold.DisplayTerms3(1,"z",Khat,G_O2,G_O3)

println("Center-Manifold Evaluation Complete.")








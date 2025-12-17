# Center Manifold Evaluation
println("Center-Manifold Evaluation on the Ginzburg Landau system.")

ifnormal = true

#if ifnormal
  include("GL_Asymptotic_CM_NormalForm2.jl")
#else
#  include("GL_Asymptotic_CM_GraphForm.jl")
#end

#CenterManifold.DisplayTerms3(1,"z",Khat,G_O2,G_O3)
CenterManifold.DisplayTerms3(1,"z",Khat,G2,G3)
CenterManifold.DisplayTerms3(2,"z",Khat,G2,G3)

println("Center-Manifold Evaluation Complete.")








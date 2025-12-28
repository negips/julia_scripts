# Setup Numerical GL

using PyPlot
#-------------------------------------------------- 

include("GL_SEM_Setup.jl")

screen = 2
Grh    = setgraphics(screen)

# GL Parameters
δ     = Set_GL_CriticalParams()
δc    = conj.(δ)

include("GL_OP_Setup.jl")

println("Ginzburg Landau Setup Done.")
















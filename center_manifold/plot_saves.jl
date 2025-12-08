# Plotting from saved files

using FileIO
using PyPlot

close("all")

dir         = "SavedFiles/"
slfname1    = "SL_resonant_Parametric.jld2"
slfname2    = "SL_nonresonant_Parametric.jld2"
slfname3    = "SL_resonant_Parametric2.jld2"
slfname4    = "SL_nonresonant_parametric2.jld2"

glfname1    = "GL_resonant_Parametric.jld2"
glfname2    = "GL_nonresonant_Parametric.jld2"
glfname3    = "GL_resonant_Parametric2.jld2"
glfname4    = "GL_nonresonant_Parametric2.jld2"

SLD1        = load(dir*slfname1) 
GLD1        = load(dir*glfname1)

GLθA1       = get(GLD1,"θA",[])
GLPA1       = get(GLD1,"Peak_Amp",[])
SLθA1       = get(SLD1,"θA",[])
SLPA1       = get(SLD1,"Peak_Amp",[])

plot(GLθA1,GLPA1,linestyle="none",marker="o")
plot(SLθA1,SLPA1,linestyle="none",marker="o",markerfacecolor="none")

# Non-Resonant
SLD2        = load(dir*slfname2) 
GLD2        = load(dir*glfname2)

GLθA2       = get(GLD2,"θA",[])
GLPA2       = get(GLD2,"Peak_Amp",[])
SLθA2       = get(SLD2,"θA",[])
SLPA2       = get(SLD2,"Peak_Amp",[])

plot(GLθA2,GLPA2,linestyle="none",marker="o")
plot(SLθA2,SLPA2,linestyle="none",marker="o",markerfacecolor="none")


println("Done.")

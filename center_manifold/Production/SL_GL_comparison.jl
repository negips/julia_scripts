# Plotting the reduced and full response curves

using FileIO
using PyPlot
using Printf

include("GL_Functions.jl")
include("../Module_CenterManifold/CenterManifold.jl")

screen = 2
Grh    = setgraphics(screen)

SLf1 = "SL_nonresonant_Parametric2.jld2"
GLf1 = "GL_nonresonant_Parametric2.jld2"

SLf2 = "SL_resonant_Parametric2.jld2"
GLf2 = "GL_resonant_Parametric2.jld2"

SLf3 = "SL_diffusion_Parametric.jld2"
GLf3 = "GL_diffusion_Parametric.jld2"

#-------------------------------------------------- 
h1    = figure(num=1,figsize=Grh.figsz3);
ax1   = gca()
ax1.cla()
ax1.set_xlabel(L"x",fontsize=Grh.lafs)
ax1.set_ylabel(L"A",fontsize=Grh.lafs)
ax1.set_title(L"Non-Resonant",fontsize=Grh.lafs)

cm    = get_cmap("tab20")

nk    = 8

# GL
#-------------------------------------------------- 
GL1    = load(GLf1)

GL1_xg = get(GL1,"xg",[])
GL1_v  = get(GL1,"vlast",[])

for ik in 1:2:nk
  ax1.plot(GL1_xg,abs.(GL1_v[:,ik]),linewidth=2,linestyle="-",color=cm(ik-1))
end  



# SL
#-------------------------------------------------- 
SL1   = load(SLf1)

SL1_xg = get(SL1,"xg",[])
Nby2   = length(SL1_xg)

SL1_θA = get(SL1,"θA",[])

SL1_M  = get(SL1,"Hist_Mode",[])
Y1     = get(SL1,"Vext",[])
Y2     = get(SL1,"Y2",[])
Y3     = get(SL1,"Y3",[])

ind1   = 1:Nby2

for ik in 1:2:nk
  lab     = @sprintf("|θ|= %.2f",SL1_θA[ik])
  z       = SL1_M[end,:,ik]
  SL1_fld = CenterManifold.GetAsymptoticField3(z,Y1,Y2,Y3)
  ax1.plot(SL1_xg,abs.(SL1_fld[ind1]),linewidth=1,linestyle="--",color=cm(ik-1),marker="o",markevery=10,label=L"%$lab")
end

ax1.legend(fontsize=Grh.lgfs)

# Resonant
#---------------------------------------------------------------------- 

nk    = 8

h2    = figure(num=2,figsize=Grh.figsz3);
ax2   = gca()
ax2.cla()
ax2.set_xlabel(L"x",fontsize=Grh.lafs)
ax2.set_ylabel(L"A",fontsize=Grh.lafs)
ax2.set_title(L"Resonant",fontsize=Grh.lafs)

# GL
#-------------------------------------------------- 
GL2    = load(GLf2)

GL2_xg = get(GL2,"xg",[])
GL2_v  = get(GL2,"vlast",[])

for ik in 1:2:nk
  ax2.plot(GL2_xg,abs.(GL2_v[:,ik]),linewidth=2,linestyle="-",color=cm(ik-1))
end  



# SL
#-------------------------------------------------- 
SL2   = load(SLf2)

SL2_xg = get(SL2,"xg",[])
SL2_θA = get(SL2,"θA",[])

SL2_M  = get(SL2,"Hist_Mode",[])
Y1Res  = get(SL2,"Vext",[])
Y2Res  = get(SL2,"Y2",[])
Y3Res  = get(SL2,"Y3",[])

for ik in 1:2:nk
  lab     = @sprintf("|θ|= %.2f",SL2_θA[ik])
  z       = SL2_M[end,:,ik]
  SL2_fld = CenterManifold.GetAsymptoticField3(z,Y1Res,Y2Res,Y3Res)
  ax2.plot(SL2_xg,abs.(SL2_fld[ind1]),linewidth=1,linestyle="--",color=cm(ik-1),marker="o",markevery=10,label=L"%$lab")
end
ax2.legend(fontsize=Grh.lgfs)



# Diffusion Coefficient Perturbation (δ4)
#---------------------------------------------------------------------- 

nk    = 4
dk    = 1

h3    = figure(num=3,figsize=Grh.figsz3);
ax3   = gca()
ax3.cla()
ax3.set_xlabel(L"x",fontsize=Grh.lafs)
ax3.set_ylabel(L"A",fontsize=Grh.lafs)
ax3.set_title(L"Diffusion",fontsize=Grh.lafs)

# GL
#-------------------------------------------------- 
GL3    = load(GLf3)

GL3_xg = get(GL3,"xg",[])
GL3_v  = get(GL3,"vlast",[])

for ik in 1:dk:nk
  ax3.plot(GL3_xg,abs.(GL3_v[:,ik]),linewidth=2,linestyle="-",color=cm(ik-1))
end  



# SL
#-------------------------------------------------- 
SL3    = load(SLf3)

SL3_xg = get(SL3,"xg",[])
SL3_δ4 = get(SL3,"δ4",[])

SL3_M  = get(SL3,"Hist_Mode",[])
Y1Dif  = get(SL3,"Vext",[])
Y2Dif  = get(SL3,"Y2",[])
Y3Dif  = get(SL3,"Y3",[])

for ik in 1:dk:nk
  lab     = @sprintf("δ4'= %.2f",SL3_δ4[ik])
  z       = SL3_M[end,:,ik]
  SL3_fld = CenterManifold.GetAsymptoticField3(z,Y1Dif,Y2Dif,Y3Dif)
  ax3.plot(SL3_xg,abs.(SL3_fld[ind1]),linewidth=1,linestyle="--",color=cm(ik-1),marker="o",markevery=10,label=L"%$lab")
end
ax3.legend(fontsize=Grh.lgfs)


println("Done.")


























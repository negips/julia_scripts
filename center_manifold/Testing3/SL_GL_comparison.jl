# Plotting the reduced and full response curves

using FileIO
using PyPlot
using Printf

include("GL_Functions.jl")
include("../Module_CenterManifold/CenterManifold.jl")

screen = 2
Grh    = setgraphics(screen)

file1 = "SL_nonresonant_Parametric1.jld2"
file2 = "SL_resonant_Parametric2.jld2"

file3 = "GL_nonresonant_Parametric1.jld2"
file4 = "GL_resonant_Parametric2.jld2"
#-------------------------------------------------- 
h1    = figure(num=1,figsize=Grh.figsz3);
ax1   = gca()
ax1.cla()
ax1.set_xlabel(L"x",fontsize=Grh.lafs)
ax1.set_ylabel(L"A",fontsize=Grh.lafs)
ax1.set_title(L"Non-Resonant",fontsize=Grh.lafs)

cm    = get_cmap("tab20")

nk    = 4

# GL
#-------------------------------------------------- 
GL1    = load(file3)
GL2    = load(file4)

GL1_xg = get(GL1,"xg",[])
GL1_v  = get(GL1,"vlast",[])

for ik in 1:nk
  ax1.plot(GL1_xg,abs.(GL1_v[:,ik]),linewidth=2,linestyle="-",color=cm(ik-1))
end  



# SL
#-------------------------------------------------- 
SL1   = load(file1)
SL2   = load(file2)

SL1_xg = get(SL1,"xg",[])
Nby2   = length(SL1_xg)

SL1_θA = get(SL1,"θA",[])

SL1_M  = get(SL1,"Hist_Mode",[])
Y1     = get(SL1,"Vext",[])
Y2     = get(SL1,"Y2",[])
Y3     = get(SL1,"Y3",[])

ind1   = 1:Nby2

for ik in 1:nk
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
GL2    = load(file4)

GL2_xg = get(GL2,"xg",[])
GL2_v  = get(GL2,"vlast",[])

for ik in 1:2:nk
  ax2.plot(GL1_xg,abs.(GL2_v[:,ik]),linewidth=2,linestyle="-",color=cm(ik-1))
end  



# SL
#-------------------------------------------------- 
SL2   = load(file2)

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


println("Done.")







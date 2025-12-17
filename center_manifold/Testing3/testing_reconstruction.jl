# Testing Field Reconstruction.

# Assuming we have Vext,Y2,Y3

h10  = figure(num=10,figsize=Grh.figsz2)
ax10 = gca()
ax10.cla()

Z     = zeros(ComplexF64,m)
phase = 1.30*π
Z[1]  = 1.0e-1*exp(phase*1.0im)
Z[2]  = Z[1]'
Z[5]  = 0.1
Z[6]  = Z[5]'

Z1    = CenterManifold.EvaluateNonLinear(Z,1)
Z2    = CenterManifold.EvaluateNonLinear(Z,2)
Z3    = CenterManifold.EvaluateNonLinear(Z,3)

Y01   = copy(Vext)
Y02   = copy(Y2)  
Y03   = copy(Y3)

fld01 = (Y01*Z1)[ind1] 
fld02 = (Y02*Z2)[ind1]
fld03 = (Y03*Z3)[ind1]

fld00 = fld01 .+ fld02 .+ fld03

ax10.plot(xg,abs.(fld01), linestyle="-", linewidth=2,color=cm(0),label="O(1)")
ax10.plot(xg,real.(fld01),linestyle="--",linewidth=1,color=cm(0),label="O(1)")
ax10.plot(xg,imag.(fld01),linestyle="-.",linewidth=1,color=cm(0),label="O(1)")

ax10.plot(xg,abs.(fld02), linestyle="-", linewidth=2,color=cm(1),label="O(2)")
ax10.plot(xg,real.(fld02),linestyle="--",linewidth=1,color=cm(1),label="O(2)")
ax10.plot(xg,imag.(fld02),linestyle="-.",linewidth=1,color=cm(1),label="O(2)")

ax10.plot(xg,abs.(fld03), linestyle="-", linewidth=2,color=cm(2),label="O(3)")
ax10.plot(xg,real.(fld03),linestyle="--",linewidth=1,color=cm(2),label="O(3)")
ax10.plot(xg,imag.(fld03),linestyle="-.",linewidth=1,color=cm(2),label="O(3)")

ax10.plot(xg,abs.(fld00), linestyle="-", linewidth=2,color=cm(3),label="Total")
ax10.plot(xg,real.(fld00),linestyle="--",linewidth=1,color=cm(3),label="Total")
ax10.plot(xg,imag.(fld00),linestyle="-.",linewidth=1,color=cm(3),label="Total")

ax10.legend(ncols=4,fontsize=Grh.lgfs)

println("Done")




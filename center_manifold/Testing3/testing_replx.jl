
include("../Module_StepperArnoldi/StepperArnoldi.jl")

λf                = λext[5]
ff                = copy(Lext[ind1,5])
ifres             = fill(false,nsys)
bb1               = [Bg; 1.0]

# Test REPLx
xt5         = copy(DArn.evecs[:,1])
dxt5        = StepperArnoldi.REPLx(xt5,OPg,Bg,VSys[ind1,:],WSys[ind2,:],σ,fill(false,nsys),ff,λf)
dxt5_B      = 0.0*dxt5
for i in 1:ndof
  dxt5_B[i] = dxt5[i]/Bg[i]
end
dxt5_B[ndof+1] = dxt5[ndof+1]/1.0
dxt5_B_xt5     = dxt5_B./xt5
#------------------------------  
xt6         = zeros(ComplexF64,ndof+1)
xt6[ndof+1] = 1.0
dxt6        = StepperArnoldi.REPLx(xt6,OPg,Bg,VSys[ind1,:],WSys[ind2,:],σ,fill(false,nsys),ff,λf)
dxt6_B      = 0.0*dxt6
for i in 1:ndof
  dxt6_B[i] = dxt6[i]/Bg[i]
end
dxt6_B[ndof+1] = dxt6[ndof+1]/1.0
dxt6_B_xt6     = dxt6_B./xt6

#------------------------------  
xt7         = zeros(ComplexF64,ndof+1)
xt7[ndof+1] = 1.0
dxt7        = zeros(ComplexF64,ndof+1)
@views StepperArnoldi.RPLx!(dxt7[1:ndof],xt7[1:ndof],OPg,Bg,VSys[ind1,:],WSys[ind2,:],σ,fill(false,nsys))
dxt7_B      = 0.0*dxt7
for i in 1:ndof
  dxt7_B[i] = dxt7[i]/Bg[i]
end
dxt7_B[ndof+1] = dxt7[ndof+1]/1.0
dxt7_B_xt7     = dxt7_B./xt7

println("Done.")

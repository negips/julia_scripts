# Testing the Laplacian calculation.

Q2  = [Q     0.0*Q;
       0.0*Q     Q]  
QT2 = [QT    0.0*QT;
       0.0*QT    QT]  

Forcing = copy(F)

h10  = figure(num=10)
ax10 = gca()
ax10.cla()

XL   = GeoM.xm1[:]
l1   = length(XL)

ind1L = 1:l1
ind2L = l1+1:2*l1

θAmp1 = θA[1]

x1    = -0.2

ax10.plot(xg,abs.(θAmp1*Vext[ind1,5]),color="blue")
ax10.plot(xg,abs.(x1*Vext[ind1,1]),color="red")

ax10.plot(xg,abs.(θAmp1*Vext[ind1,5] + x1*Vext[ind1,1]),color="black")

vtmp1 = 0.0*ψ
θtmp1 = [vt(θAmp1)]

df,dθ    = FNGL(vtmp1,θtmp1)
#ax10.plot(xg,real.(df))

OP2_RK4!(FNGL,vtmp1,θtmp1,dt,vwork,θwork)

#ax10.plot(xg,real.(vtmp1))

iters = 100000

for i in 1:iters
  global vtmp1,θtmp1

  OP2_RK4!(FNGL,vtmp1,θtmp1,dt,vwork,θwork)
end  

ax10.plot(xg,real.(vtmp1),linestyle="--")
ax10.plot(xg,abs.(vtmp1),linestyle="-",linewidth=2)

println("Done")







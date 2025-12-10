# Testing the Laplacian calculation.


ϕ1 = Vext[:,1]
ϕ5 = Vext[:,5]

Q2  = [Q     0.0*Q;
       0.0*Q     Q]  
QT2 = [QT    0.0*QT;
       0.0*QT    QT]  

ϕ1L = [Q*ϕ1[ind1]; Q*ϕ1[ind2]]
ϕ5L = [Q*ϕ5[ind1]; Q*ϕ5[ind2]]

h10  = figure(num=10)
ax10 = gca()
ax10.cla()

XL   = GeoM.xm1[:]
l1   = length(XL)

ind1L = 1:l1
ind2L = l1+1:2*l1

ax10.plot(XL,real.(ϕ1L[ind1L]))
ax10.plot(XL,real.(ϕ5L[ind1L]))

Lapϕ1L = Lap*ϕ1L[ind1L]
Lapϕ5L = Lap*ϕ5L[ind1L]

ax10.plot(XL,real.(Lapϕ1L))
ax10.plot(XL,real.(Lapϕ5L))

@views SEM1D.SEM_SetBC!(Lapϕ1L,Inp.lbc,Inp.rbc)
@views SEM1D.SEM_SetBC!(Lapϕ5L,Inp.lbc,Inp.rbc)

BiLapϕ1L = Q*(Bgi.*(QT*Lapϕ1L))
BiLapϕ5L = Q*(Bgi.*(QT*Lapϕ5L))

ax10.plot(XL,real.(BiLapϕ1L))
ax10.plot(XL,real.(BiLapϕ5L))

println("Done")




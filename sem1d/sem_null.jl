println("Eigen spectrum of convection diffusion operator")

using LinearAlgebra
using PyPlot,PyCall

close("all")

# Include the function files
include("sem_main.jl")

r,c = size(L);
Lsub = L[2:r,2:c];

F  = eigen(Lsub);

F.values

lr = real.(F.values);
li = imag.(F.values);

ii = sortperm(lr,rev=true);

figure(num=1);
plot(li,lr,linestyle="none",marker="o",markersize=8)

xall  = zeros(Float64,r);
for i in 1:nel
   global xall   
   j1 = (i-1)*N + 1;
   j2 = j1+lx1-1;
   xall[j1:j2] = Geom.xm1[:,i];
end   


#figure(num=2)
#plot(xall[2:end],real.(F.vectors[:,ii[1]]));
#plot(xall[2:end],imag.(F.vectors[:,ii[1]]));

r,c = size(Lsub);

eye = zeros(Complex,r,r);
for i in 1:r
  eye[i,i] = 1.;
end  

A = Complex.(Lsub);

i = sqrt(Complex(-1.));

#ω = -1.0 + 0.1*i*ω;

figure(num=2)
for j in (3) #30
  global pl, ns
  
  ω = F.values[ii[j]];

  println(ω)

  D = sqrt(A);

  ω1 = sqrt(ω); 
  B  = D + i*ω1.*eye;
  ns = nullspace(B)

  F2 = eigen(D);
  ns = F2.vectors;
 
  r2,c2=size(ns);
  println("$c2 null vectors found")

  figure(num=1);
  lr2 = real.(F2.values);
  li2 = imag.(F2.values);
  
  ii2 = sortperm(lr2,rev=true);
  
  figure(num=1);
  plot(li2,lr2,linestyle="none",marker=".",markersize=8)

#  for k in 1:c2
#    if k>1
#      pl[1].remove();
#    end  
#       
#    pl = plot(xall[2:end],real.(ns[:,k]),marker=".");
#    pause(2.)
#  end  
end  
#plot(xall[2:end],imag.(F.vectors[:,ii[1]]));




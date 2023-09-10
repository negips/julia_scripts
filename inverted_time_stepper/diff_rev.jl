
# Time Reversal
println("Time Reversal")
plotupd = 1;

deltastep = 5

for i in nsteps-1:-1:(nsteps-deltastep)
  global V,Vlag,Rlag
  global t
  global pl

  t = t - dt;

  bdf = bdf3;
  ext = ex2;
     
  if i==1
    bdf = bdf1;
    ext = ex0;
  elseif i==2  
    bdf = bdf2;
    ext = ex1;
#  elseif i==3 
#    bdf = bdf3;
#    ext = ex2;
  end

  A = Lg;
#  A[1,:] = zeros(Float64,1,r);
#  A[1,1] = 1.

# No Extrapolate A*x
  rhs     = 0.0*V;
  rhs1    = ext[2]*rhs + ext[3]*Rlag[:,1] + ext[4]*Rlag[:,2];

  Rlag[:,2] = Rlag[:,1];
  Rlag[:,1] = rhs;

  bdlag    = 1. /dt*(bdf[1]*V + bdf[2]*Vlag[:,1] + bdf[3]*Vlag[:,2]) # + bdf[4]*Vlag[:,2]; 
  Rhs      = rhs1 + bdlag - A*V;

  

#  M        = bdf[4]/dt*I # + A
#  S        = -gmres(M,Rhs,abstol=1.0e-10)
   S       = -Rhs*dt/bdf[4]   

#  Vlag[:,3] = Vlag[:,2];
#  Vlag[:,2] = Vlag[:,1];
#  Vlag[:,1] = V;

  V          = Vlag[:,1]
  Vlag[:,1]  = Vlag[:,2]
  Vlag[:,2]  = S

  if mod(i,plotupd)==0
    if (i>plotupd)
       pl[1].remove();
    end   
  
    pl = plot(xall,Vlag[:,2],color=rgba1);
 
    pause(0.001)
  end  

end














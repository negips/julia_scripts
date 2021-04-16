println("Standing wave as a result of forward and backward traveling waves")

using PyPlot,PyCall

close("all")

nsteps      = 10000;
plotstep    = 10;
dt          = 0.01;

x = range(0,stop=2*pi,length=1000);
i     = sqrt(Complex(-1));

sigma1 = 0.0;
omega1 = 4.0;
ki1    = 4.0;
kr1    = 0.0;

sigma2 = 0.0;
omega2 = 1.0;
ki2    = 0.5;
kr2    = 0.0;


cm    = get_cmap("tab10");
rgba0 = cm(0); 
rgba1 = cm(1);
rgba2 = cm(2);

ϕ     = 0; # 0*rand(Float64);
a1    = 0.1;
a2    = 1.0;

for j in 1:nsteps
   global plf,plb,pls

   if (j>plotstep) & (mod(j,plotstep)==0)
#     plf[1].remove();
#     plb[1].remove();
     pls[1].remove();
   end  

   t  = j*dt; 
   sw11 = a1*sin.(ki1.*x).*sin.(omega1*t );

   sw21 = a2*sin.(ki2.*x).*sin.(omega2*t  .+ ϕ);

   sw1r = real.(sw11);
   sw2r = real.(sw21);

   trw  = sw1r + sw2r;

   if mod(j,plotstep)==0
#     plf = plot(x,sw1r,linestyle="--",color=rgba0)
#     plb = plot(x,sw2r,linestyle="--",color=rgba1)
     pls = plot(x,trw,linestyle="-",color=rgba2)
     ax  = gca();
#     ax.set_ylim([-3.5,3.5]);
     
     pause(0.001)
   end  

end   




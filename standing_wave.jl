println("Standing wave as a result of forward and backward traveling waves")

using PyPlot,PyCall

close("all")

nsteps      = 10000;
plotstep    = 10;
dt          = 0.01;

x = range(0,stop=20*pi,length=10000);

i     = sqrt(Complex(-1));
sigma = 0.01;
omega = 1.0;
ki    = 1.0;
kr    = 0.05;

cm    = get_cmap("tab10");
rgba0 = cm(0); 
rgba1 = cm(1);
rgba2 = cm(2);

ϕ     = 0*rand(Float64);
a1    = 1.0;
a2    = 1.0;

for j in 1:nsteps
   global plf,plb,pls

   if (j>plotstep) & (mod(j,plotstep)==0)
     plf[1].remove();
     plb[1].remove();
     pls[1].remove();
   end  

   t  = j*dt; 
   fw = a1*exp.((kr + i*ki).*x .+ sigma*t .- i*omega*t );
   bw = a2*exp.((kr + i*ki).*x .+ sigma*t .+ i*omega*t  .+ ϕ);

   fwr = real.(fw);
   bwr = real.(bw);

   swr = fwr + bwr;

   if mod(j,plotstep)==0
     plf = plot(x,fwr,linestyle="--",color=rgba0)
     plb = plot(x,bwr,linestyle="--",color=rgba1)
     pls = plot(x,swr,linestyle="-",color=rgba2)
     ax  = gca();
#     ax.set_ylim([-3.5,3.5]);
     
     pause(0.001)
   end  

end   




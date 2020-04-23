#!/usr/bin/julia
# Source code for Homework 3

println("Homework 2, Exercise 2c")

#import Pkg
#Pkg.add("PyPlot")
#Pkg.add("PyCall")
#Pkg.add("Blink")

using Blink,Random
using PyPlot,PyCall,LinearAlgebra,SparseArrays
using LaTeXStrings,TimerOutputs

include("GS_npass.jl")
include("GS_modified.jl")
include("mygmres.jl")

const TO = TimerOutput();

close("all")      # close all open figures

ifmod = false;    # if we use modified GS
ngs = 1;          # No of GS passes

errtol = 1.0e-10; # error tolerance for sub diagonal elements

cols = ["b","r","k","g"];
lafs=20;
lgfs=16;

m = 500;
#A = spzeros(m,m);

alpha_range=[1,100];
nalpha=length(alpha_range);

ifploterr=false;
ifploteig=false;

for ii = 1:nalpha

  alpha = alpha_range[ii];
  
  Random.seed!(5);
  A = sprand(Float64,m,m,0.5);
  A = A + alpha*sparse(1.0I,m,m);
  A = A/norm(A);
  b = rand(Float64,m,1);
  
  xexact = @timeit TO "Backslash.$alpha" A\b;
  resid = norm(A*xexact - b);
  println("Exact resid=$resid for alpha=$alpha")
  
  tol = 1.0e-14;
  maxits = 100;
  ifcomplex = false
# Dummy so I get timings right
  x1,Q1,H1,r1,xerr1 = mygmres(A,b,ngs,ifmod,tol,maxits,xexact,ifcomplex); 
  x,Q,H,r,xerr = @timeit TO "MyGMRES.$alpha" mygmres(A,b,ngs,ifmod,tol,maxits,xexact,ifcomplex);

  resid=r[end];
#  println(size(Q))
  println("GMRES resid=$resid for alpha=$alpha for m=$maxits iterations")

  str1=latexstring("\$||x - x_{*}||; \\alpha=$alpha\$")
  str2=latexstring("\$||Ax - b||; \\alpha=$alpha\$")
 
  if (ifploterr)
    pl1 = semilogy(xerr,color=cols[ii],label=str1);
    pl2 = semilogy(r,color=cols[ii],linestyle="--", label=str2);
    ylabel(L"\varepsilon",fontsize=lafs)
    xlabel("Iterations",fontsize=lafs)
    fm1 = plt.get_current_fig_manager();
    fm1.window.wm_geometry("1200x800")
  end  
  
  #plt[:figure](1)
  #fm1 = plt[:get_current_fig_manager]();
  #fm1[:window][:wm_geometry]("800x600")
  
 
  evals = eigvals(Matrix(A));
  
  er = real(evals);
  ei = imag(evals);
  
  rmax,imax = findmax(er);
  rmin,imin = findmin(er);
  
  cen = (rmax+rmin)/2;
  rad = rmax-cen;
  cfac = rad/cen;
  theta = collect(0.:0.00001:2*pi);
  x1=cen .+ rad*cos.(theta);
  y1=0. .+ rad*sin.(theta);
  
  niters = length(xerr);
  ik = collect(1:niters);
  conv1 = 10. *r[1]*(cfac).^ik;
  #pl3 = semilogy(conv1,color=cols[1],linestyle=":",label="Theoretical")
  
  
  # Creating a more sophisticated circle
  ind = sortperm(er);
  ers = er[ind];
  eis = ei[ind];
  neig = length(ers);
  rmax = ers[neig-1];
  rmin = ers[1];
  l1 = rmax - rmin;
  l2 = abs(eis[2]);
  tanb = l2/l1;
  l3 = l1*tanb^2;
  dia = l1 + l3;
  opp_x = rmax - dia;
  cen = (rmax + opp_x)/2;
  rad = dia/2;
  
  c2 = ers[neig];
  
  rho = sqrt(rad*(rad + abs(cen-c2))/abs(cen*c2))
  
  cfac = rad/cen;
  theta = collect(0.:0.00001:2*pi);
  x2=cen .+ rad*cos.(theta);
  y2=0. .+ rad*sin.(theta);
  conv2 = 10000. *r[1]*(cfac).^(ik);
  conv3 = 1000. *r[1]*(rho).^(ik.-1);

  if (ifploterr)
    if alpha!=1
      str3=latexstring("\$(r/c)^{m}; \\alpha=$alpha\$")       
      pl3 = semilogy(conv2,color=cols[ii],linestyle=":",label=str3)
     #pl4 = semilogy(conv3,color=cols[1],linestyle="-.",label=L"\rho^{m}")
    end
    xlim(0,150)
  end  
  

  if (ifploteig)
    #fm2 = plt[:get_current_fig_manager]();
    fm2 = plt.get_current_fig_manager();
    #fm2[:window][:wm_geometry]("800x600")
    fm2.window.wm_geometry("800x600")
    str4=latexstring("\$\\alpha=$alpha\$")       
    plot(er,ei,linestyle="none",marker=".",markersize=4,label=str4)
    plot(x2,y2,linewidth=2)
    ax=gca();
    #ax[:set_position]([0.15, 0.15, 0.8, 0.7])
    ax.set_position([0.15, 0.15, 0.8, 0.7])
    ylabel(L"\lambda_{i}",fontsize=lafs)
    xlabel(L"\lambda_{r}",fontsize=lafs)
    #tight_layout
  end

end

if (ifploterr)
  legend(loc="best",fontsize=lgfs)
  fname="ex2a_err$alpha.png"
  savefig(fname)
end

if (ifploteig)
  legend(loc="best",fontsize=lgfs)
  fname="ex2a_eigs$alpha.png"
  savefig(fname)
end  

TO

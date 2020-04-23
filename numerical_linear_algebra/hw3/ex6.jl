using LinearAlgebra, Random
using Plots

gr()
#epsi = 0.5;
#A = [pi 1;
#     0 pi+epsi];
#
#X  = [1 1/sqrt(1+epsi^2);
#     0 sqrt(epsi^2/(1+epsi^2))];
#D  = [pi 0;
#     0 pi+epsi];
#Xi = [1 -1/epsi;
#      0 sqrt((1+epsi^2)/(epsi^2))];

# The exact expression for the exp(A)

epsiv = exp10.(-2:-0.1:-10);
#epsiv = 1;
#expA(epsiVal)=[exp(pi) (exp(pi+epsiVal)-exp(pi))/epsiVal;
#            0 exp(pi+epsiVal)];

fexp = (z)->exp(z);
errorv = Array{Float64}(undef,0);

for epsi in epsiv

    A = [pi 1;
         0 pi+epsi];
    
    D,V = eigen(A);
    F   = diagm(0=>exp.(D));
    FJ  = V*F*inv(V);
    
    #    FE = expA(epsi);
    alfa = fexp(pi)-(fexp(pi+epsi)-fexp(pi) )*pi/epsi;
    beta = (fexp(pi+epsi)-fexp(pi))/epsi;
    FE = [alfa+beta*pi beta;
          0 alfa+beta*(pi+epsi)];
    
    errJ = norm(FE-FJ);
    append!(errorv,errJ);
end

p1 = Plots.plot(epsiv,errorv,marker=2,xaxis=:log10,yaxis=:log10,label="")
xlabel!("\\epsilon"); ylabel!("||f(A)-F||")

fname = "Ex6_epsi_err.pdf";
Plots.savefig(p1,fname)
display(p1)

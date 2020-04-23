# Load the necessary modules
using LinearAlgebra
using Plots
using MatrixDepot, Random, Arpack
using MAT, LaTeXStrings


########################################
#           EXERCISE 2                 #
########################################

include("functionsEx2.jl")

# Matrix from HW1 exercise 2
A    = [1 2 3;2 2 2; 3 2 9];

x0   = [1;1;1];     # Initial guess
maxk = 20;          # Max iterations

Avals, Avecs = eigen(A);
maxAvals = maximum(Avals);
sort!(Avals);
p99 = scatter(Avals,marker=2)

maxEigPM,vPM,itEigPM=powerIt(A,x0,maxk);
itErrPM = norm.(itEigPM.-maxAvals,2);

# There is a problem with the error being exactly 0 at some points
# presenting an issue with the log plot. The next lines are done
# to fix the error to numerical epsilon if the eigErr is too low

for j=1:length(itErrPM)
    itErrPM[j]=max(itErrPM[j],eps())
end

gr()
p1=plot(collect(1:1:maxk),itErrPM,label="eigErr in PM",yaxis=:log10,marker=2)
#p1=plot!(p1,collect(1:1:maxk),itErrRQI,label="eigErr in RQI",yaxis=:log10,marker=2)

xlabel!("Iterations")
ylabel!("\$ abs(\\lambda-\\lambda^k) \$")
ylims!((1e-16,1e1))
savefig(p1,"errLambda_PM.pdf")

maxk = 10;          # Max iterations

maxEigRQI,vRQI,itEigRQI=RQit(A,x0,maxk);
itErrRQI = norm.(itEigRQI.-maxAvals,2);

# There is a problem with the error being exactly 0 at some points
# presenting an issue with the log plot. The next lines are done
# to fix the error to numerical epsilon if the eigErr is too low

for j=1:length(itErrRQI)
    itErrRQI[j]=max(itErrRQI[j],eps())
end

p2=plot(collect(1:1:maxk),itErrRQI,label="eigErr in RQI",yaxis=:log10,marker=2)

xlabel!("Iterations")
ylabel!("\$ \\lambda-\\lambda^k \$")
ylims!((1e-16,1e1))
savefig(p2,"errLambda_RQI.pdf")

# Matrix from HW1 exercise 2
A    = [1 2 4;2 2 2; 3 2 9];
x0   = [1;1;1];     # Initial guess

Avals, Avecs = eigen(A);
maxAvals = maximum(Avals);
sort!(Avals);

maxEigRQI2,vRQI2,itEigRQI2=RQit(A,x0,maxk);
itErrRQI2 = norm.(itEigRQI2.-maxAvals,2);

# There is a problem with the error being exactly 0 at some points
# presenting an issue with the log plot. The next lines are done
# to fix the error to numerical epsilon if the eigErr is too low

for j=1:length(itErrRQI2)
    itErrRQI2[j]=max(itErrRQI2[j],eps())
end

p3=plot(collect(1:1:maxk),itErrRQI,yaxis=:log10,marker=2)
p3=plot!(p3,collect(1:1:maxk),itErrRQI2,yaxis=:log10,marker=2)
labels = ["\$ A= A^T \$" "\$ A\\neq A^T \$"];
p3 = plot(p3,label=labels);
display(p3)

xlabel!("Iterations")
ylabel!("\$ abs(\\lambda-\\lambda^k) \$")
ylims!((1e-16,1e1))
savefig(p3,"errLambda_RQI_comp.pdf")


#gr()
#p2=plot(collect(1:1:maxk),itErrPM,label="eigErr in PM",yaxis=:log10,marker=2)
#p2=plot!(p2,collect(1:1:maxk),itErrRQI,label="eigErr in RQI",yaxis=:log10,marker=2)

#xlabel!("Iterations")
#ylabel!("\$ \\lambda-\\lambda^k \$")
#ylims!((1e-16,1e1))

using LinearAlgebra
using Plots
using MatrixDepot, Random, Arpack
using MAT, Printf, Statistics

println("Exercise 6")
include("arnoldi.jl")

gr()
B      = matread("Bwedge.mat")["B"];
λB     = matread("Bwedge.mat")["B_eigvals"];
isortB = sortperm(λB[:,1].-mean(λB[:,1]),by=abs);
λB     = λB[isortB];

p4 = scatter(real(λB),imag(λB),marker=2,color=:red)

C   = λB[end];
rho = 0.5; 

# Plot the circle
theta   = collect(0:0.01:2*pi);
xcircle = rho*cos.(theta).+real(C);
ycircle = rho*sin.(theta).+imag(C);

p4 = plot!(p4,xcircle,ycircle,color=:red)

# Compute now the convergence and the bound for the convergence ratio

C   = 0.5*(λB[end-1]+λB[end-2]);
rho = maximum(abs.(λB[λB .!= λB[end]].-C) ); 

# Plot the circle
theta   = collect(0:0.01:2*pi);
xcircle = rho*cos.(theta).+real(C);
ycircle = rho*sin.(theta).+imag(C);

p4 = plot!(p4,xcircle,ycircle,line=(:dash,2),color=:blue)
p4 = scatter!(p4,[real(C)],[imag(C)],marker=(:+,4),color=:blue)
xlabel!("Real")
ylabel!("Imag")

labels = ["Eigs of B" "Favoured eig" "\$ \\rho = 19.38 , C = 1.15+0.38im, cfact = 0.402 \$" ""]
display(plot(p4,label=labels,legend=:topleft))
savefig(p4,"exercise6_ab.pdf")

cfactor = abs(rho/(λB[end]-C));

# Do in a loop to generate the figure
 mvec  = [2; 4; 8; 10; 20; 30; 40];
# mvec = [2];
λmax  = zeros(Complex{Float64},length(mvec));
ebounds = zeros(length(mvec));

p = 0;

# Restart the random seed
Random.seed!(0);
b = rand(length(B[:,1]),1);
#b = B[:,1];

# Define if there will be a shift
#mu = -11+2im;
mu = 0;
#B = inv(B-mu*I);



for m in mvec
    println("m ",m)
    global p += 1;
    
    # Do the Arnoldi
    Q,H= arnoldi(B,b,m);
    λH, ϕH  = eigen(H[1:m,1:m]);
    isort   = sortperm(λH,by=abs);
    λH =  λH[isort];

    # Reverse the eigenvalue transformation if shift
    #λH = 1 ./λH.+mu;
    
    p6 = scatter!(p4,real(λH),imag(λH),marker=(:hexagon),color=:green,reuse=false)
    xlabel!("Real")
    ylabel!("Imag")
    
    labels = ["Eigs of B" "Favoured eig" "\$ \\rho = 19.38 , C = 1.15+0.38im, cfact = 0.402 \$" "" "Eigs of Hm" "" "" "" "" "" "" "" "" "" "" "" "" "" "" ""];
    p6 = plot(p6,label=labels,legend=:topleft,reuse=false);
    
    fname = @sprintf("%s%04d.pdf","exercise6_",m);
    savefig(p6,fname)

    
#    display(plot(p6))
#    display(plot(p4))
#    readline(stdin)

    λmax[p] = λH[end];
    ebounds[p] = cfactor^(m-1);
end

  
labels = ["Eigs" "" "" "" "" "Eig Arnoldi"] 

# The error in Arnoldi is normalized so that the eigenvalue is norm to 1

earnoldi = abs.( (λmax .-λB[end])./λB[end]);

p5 = plot(mvec,ebounds,yaxis=:log10,marker=2)
p5 = plot!(mvec,earnoldi,yaxis=:log10,ylims=(1e-15,1),marker=2);
xlabel!("Iterations")
ylabel!("\$ abs(\\lambda-\\lambda^k) \$")
labels = ["Error bound" "Error eig Arnoldi"]
display(plot(p5,label=labels))
savefig(p5,"exercise6_convergence.pdf")


######################################################
# Do the same for the inverted matrix

include("arnoldi_shift.jl")

gr()
B      = matread("Bwedge.mat")["B"];
mu     = -11+2im;
B2     = inv(B-mu*I);
λB2     = eigvals(B2);
#λB     = matread("Bwedge.mat")["B_eigvals"];
isortB = sortperm(λB2.-mean(λB2),by=abs);
λB2    = λB2[isortB];

p4 = scatter(real(λB2),imag(λB2),marker=2,color=:red)

C   = λB2[end];
rho = 0.05; 

# Plot the circle
theta   = collect(0:0.01:2*pi);
xcircle = rho*cos.(theta).+real(C);
ycircle = rho*sin.(theta).+imag(C);

p4 = plot!(p4,xcircle,ycircle,color=:red)
# Compute now the convergence and the bound for the convergence ratio

C2   = 0.5*(λB2[end-2]+λB2[end-3]);
rho2 = maximum(abs.(λB2[λB2 .!= λB2[end]] .- C2) ); 
rho2 = rho2+0.05*rho2;
# Plot the circle
theta   = collect(0:0.01:2*pi);
xcircle = rho2*cos.(theta).+real(C2);
ycircle = rho2*sin.(theta).+imag(C2);

p4 = plot!(p4,xcircle,ycircle,line=(:dash,2),color=:blue)
p4 = scatter!(p4,[real(C2)],[imag(C2)],marker=(:+,4),color=:blue)
xrang = xlims(p4);
yrang = ylims(p4);
xlabel!("Real")
ylabel!("Imag")
xlims!(xrang)
ylims!(yrang)

labels = ["Eigs of B" "Favoured eig" "\$ \\rho = 1.89 , C = -0.642+0.533im, cfact = 0.352 \$" ""]
display(plot(p4,label=labels,legend=:bottomright))
savefig(p4,"exercise6_d.pdf")

cfactor2 = rho2/abs((λB2[end]-C2));

# Do in a loop to generate the figure
 mvec  = [2; 4; 8; 10; 20; 30; 40];
# mvec = [2];
λmax  = zeros(Complex{Float64},length(mvec));
ebounds2 = zeros(length(mvec));
earnoldi2 = zeros(length(mvec));

p = 0;

# Restart the random seed
Random.seed!(0);
b = rand(length(B[:,1]),1)+rand(length(B[:,1]),1)im;
#b = B2[:,1];

# Define if there will be a shift
mu = -11+2im;
#mu = 0;
#B = inv(B-mu*I);
p=0;
for m in mvec
    println("m ",m)
    global p += 1;
    
    # Do the Arnoldi
    Q2,H2= arnoldi_shift(B,b,mu,m);
#    Q2,H2= arnoldi(B,b,m);
    λH2, ϕH2  = eigen(H2[1:m,1:m]);
#    isort2   = sortperm(λH2,by=abs);
#    λH2 =  λH2[isort2];

    # Reverse the eigenvalue transformation if shift
    λH2rev = 1 ./λH2.+mu;
    
    p6 = scatter!(p4,real(λH2),imag(λH2),marker=(:hexagon),color=:green,reuse=false)
    xlabel!("Real")
    ylabel!("Imag")
    xlims!(xrang)
    ylims!(yrang)
    labels = ["Eigs of B" "Favoured eig" "\$  \\rho = 1.89 , C = -0.642+0.533im, cfact = 0.352 \$" "" "Eigs of Hm" "" "" "" "" "" "" "" "" "" "" "" "" "" "" ""];
    p6 = plot(p6,label=labels,legend=:bottomright,reuse=false);    
    fname = @sprintf("%s%04d.pdf","exercise6_shift_",m);
    savefig(p6,fname)

    
#    display(plot(p6))
#    display(plot(p4))
#    readline(stdin)

    earnoldi2[p] = minimum(abs.(λH2 .-λB2[end]))./abs(λB2[end]);
    ebounds2[p] = cfactor2^(m-1);
end

  
labels = ["Eigs" "" "" "" "" "Eig Arnoldi"] 

# The error in Arnoldi is normalized so that the eigenvalue is norm to 1

#earnoldi = abs.( (λmax .-λB[end])./λB[end]);

p5 = plot(mvec,ebounds2,yaxis=:log10,marker=2)
p5 = plot!(mvec,earnoldi2,yaxis=:log10,ylims=(1e-15,1),marker=2);
xlabel!("Iterations")
ylabel!("\$ abs(\\lambda-\\lambda^k) \$")
labels = ["Error bound" "Error eig Arnoldi"]
display(plot(p5,label=labels))
savefig(p5,"exercise6_shift_convergence.pdf")

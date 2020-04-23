"""

Example for converging to the desired eigenvalue even when they are not distinct in modulus
using PM, based on the inverse iterations. The spiral matrix will be used as an example

"""

# Load the necessary modules
using LinearAlgebra
using Plots
using MatrixDepot, Random, Arpack
using MAT

function IPM(A,x0,mu,maxk)

    # Input
    # A   :  Matrix on which the RQit is used
    # x0  :  Initial guess for the eigvec
    # mu  :  Shift for inverse iterations
    # maxk:  Maximun number of iterations

    maxEig = nothing;                    # Initializes the output
    itEig  = Array{Complex{Float64}}(undef,0);    # Initializes the output

    v = x0; # Initialize the eigvec
    for k in 1:maxk
        w  = (A-mu*I)\v;
        v  = w/norm(w,2);
        maxEig = v'A*v/(v'v);
        append!(itEig,maxEig);
    end

    return maxEig,v,itEig

end

function load_spiral_matrix(n)

    # Input
    # n    : Numer of points

    tv=LinRange(0.0,2*pi,n);
    zv=zv=3 .+((1 ./(1.5 .-tv ./(2*pi))).^5).*(sin.(5 .*tv).+0.5im*cos.(5 .*tv)) ./30;
    m=length(zv);
    Random.seed!(0);
    A=randn(m,m);
    A=A\diagm(0 => zv)*A;

    return A
end

n = 100;
A = load_spiral_matrix(n);
eigval,eigvec = eigen(A);

gr()
p1 = scatter(real(eigval),imag(eigval),marker=2,xlabel="real",ylabel="imag");
display(plot(p1,legend=false))

Random.seed!(0);
epsi = 1e-12;
x0 = epsi.*rand(n,1)+rand(n,1)im;
mu = 3.0+0.38im;
maxk = 250;

maxEig,v,itEig = IPM(A,x0,mu,maxk);

p1 = scatter!(real(itEig),imag(itEig),marker=2,xlabel="real",ylabel="imag")
display(plot(p1))

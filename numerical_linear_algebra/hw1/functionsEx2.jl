"""
Functions to run on the exercise 2 HW1
Miguel Beneitez - 05112018
"""

function powerIt(A,v,maxk)

    # Input
    # A   :  Matrix on which the PM is used
    # maxk:  Maximun number of iterations
    # x0  :  Initial guess for the eigvec

    maxEig = nothing                    # Initializes the output
    itEig  = Array{Float64}(undef,0)    # Initializes the output
    
    # Power method implementation

    for k in 1:maxk
        w = A*v
        v = w/norm(w,2)
        maxEig = (v'A*v)/(v'v)
        append!(itEig,maxEig)
    end
    
    return maxEig,v,itEig

end

function RQit(A,v,maxk)

    # Input
    # A   :  Matrix on which the RQit is used
    # maxk:  Maximun number of iterations
    # x0  :  Initial guess for the eigvec

    maxEig = nothing                    # Initializes the output
    itEig  = Array{Float64}(undef,0)    # Initializes the output
    
    # Rayleigh quotient iteration
    # mu  : Internal variable corresponding to the Rayleigh quotient
    
    mu = v'A*v/norm(v'*v,2) 
    for k in 1:maxk
        w  = (A-mu*I)\v
        v  = w/norm(w,2)
        maxEig = v'A*v/(v'v)
        append!(itEig,maxEig)
        mu = maxEig
    end

    return maxEig,v,itEig

end

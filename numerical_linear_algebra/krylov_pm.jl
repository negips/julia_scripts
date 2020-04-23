function krylov_pm(A,b,m)
#   Generate Krylov space using the Power method

    n=length(b);
    Km=zeros(n,m);
    Km[:,1]=b/norm(b);

    Ab = b;  

    for k=2:m
        Ab = A*(Ab);
        Km[:,k]=Ab./(norm(Ab)); 
    end

    return Km
end

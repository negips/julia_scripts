function schur_parlett(A,f)
# The Schur-Parlett method. The function handle f is a scalar
# valued function

#    global T
#    global F  

    T,Q=schur(Matrix{ComplexF64}(A));   # complex schur form since to make
                                        # it work for complex complex eigenvalues
    n=size(A,1);
    F=zeros(ComplexF64,n,n)
    for i=1:n
        F[i,i]=f(T[i,i]);
    end
    for p=1:n-1
        for i=1:n-p
            j=i+p
#            global s
            s=T[i,j]*(F[j,j]-F[i,i])

            for k=i+1:j-1
                  s=s+T[i,k]*F[k,j]-F[i,k]*T[k,j];
            end
            F[i,j]=s/(T[j,j]-T[i,i]);
        end    
    end
    F=Q*F*Q';
    return F
end

#A=rand(3,3);
#A=A;           # Testing. Only real eigenvalues.
#f=z->sin(z)
#F=schur_parlett(A,f)

function chebdif(N,M)

#   using LinearAlgebra,ToeplitzMatrices

#  Inputs
#  N - No of Chebyshev points
#  M - no of Derivatives

#  Outputs
#  x - chebyshev points
#  DM - Chebyshev differentiation matrix


   N1 = N-1;
   eye=Matrix{Float64}(I, N, N);
   n1 = trunc(Int, floor(N/2));
   n2 = trunc(Int, ceil(N/2));
   k  = collect(0:N1);
   th = k.*pi/(N1);

   l = collect(N1:-2:-N1);
   x = sin.(pi.*l./(2*N1));

   T  = repeat(th./2,1,N);
   DX = 2*sin.(T'+T).*sin.(T'-T);
   DX = [DX[1:n1,:]; -reverse(reverse(DX[1:n2,:],dims=2),dims=1)];
   for i in 1:N
     DX[i,i]=1.
   end  

   C = real.(Toeplitz((-1.).^k,(-1.).^k));
   C[1,:] = C[1,:].*2;
   C[N,:] = C[N,:].*2;
   C[:,1] = C[:,1]./2;
   C[:,N] = C[:,N]./2;

   Z = 1 ./DX;
   for i in 1:N
     Z[i,i]=0.
   end

   D  = Matrix{Float64}(I,N,N);
   DM = Array{Float64}(undef,N,N,M);
  
   for ell in 1:M
     d = diag(D);
     D = ell.*Z.*(C.*repeat(d,1,N) - D);
     for j in 1:N
       D[j,j] = -sum(D[j,:]);
     end
     DM[:,:,ell] = D;
   end   

   return x,DM
end




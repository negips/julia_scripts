function GS_modified(Q,w,k,ifcomplex)

#  Inputs
#  Q - Matrix of orthogonal k orthogonal vectors
#  k - no of orthogonal vectors
#  w - vector to orthogonalize

#  Outputs
#  h - components along Qk vectors
#  beta - component orthogonal to Qk
#  wortho - vector orthogonal to Qk

   rw,cl=size(Q);
   if (ifcomplex)
     h = zeros(Complex,k,1);  # Initialize to 0
   else         
     h = zeros(Float64,k,1);  # Initialize to 0
   end  
      
   for i in 1:k
     g = Q[:,i]'*w;
     w = w .- Q[:,i].*g;
     h[i] = g; 
   end   

   wnorm = norm(w);
   qnew = w./wnorm;
   beta = w'*qnew;
   wortho = beta.*qnew;
      
#  Check orthogonality of qnew
#   test = norm(Q[:,1:k]'*qnew);
#   println("Orthogonality of new Vector: $test")

   return h,beta,wortho
end


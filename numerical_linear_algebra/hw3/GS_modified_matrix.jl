function GS_modified_matrix(A,k)

#  Inputs
#  A - Matrix to orthogonalize
#  k - no of columns to orthogonalize

#  Outputs
#  R - R matrix
#  Q - Orthogonal matrix

   rw,cl=size(A);
   Q = copy(A[:,1:k]);
   v = zeros(Float64,rw,1);
   R = zeros(Float64,k,k); 

   for i in 1:k

     v = Q[:,i]; 
     R[i,i] = norm(v);
     v = v./R[i,i];
     Q[:,i] = v;        # Maybe this operation is not necessary

     for j in i+1:k
       R[i,j] = v'*Q[:,j]; 
       Q[:,j] = Q[:,j] - R[i,j].*v;
     end 
   end   
      
#  Check orthogonality of Q
#   test = norm(Q[:,1:k]'*Q[:,1:k] - I);
#   println("Orthogonality of new Vector: $test")

   return Q,R
end


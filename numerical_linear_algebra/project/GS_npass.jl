function GS_npass(Q,w,k,ngs)
#  write your function

#  Inputs
#  Q - Matrix of orthonormal vectors
#  k - no of orthonormal vectors
#  w - vector to orthogonalize
#  ngs - No of Gram-Schmidt loops

#  Outputs
#  h - components along Qk vectors
#  beta - component orthogonal to Qk
#  wortho - vector orthogonal to Qk

   rw,cl=size(Q);
   h = fill(0.,k);  # Initialize to 0

   for igs in 1:ngs
#     global h
#     global w 

     g = Q[:,1:k]'*w;
     h = h .+ g; 
     w = w .- Q[:,1:k]*g;
   end   

   wnorm = norm(w);
   qnew = w./wnorm;
   beta = wnorm;
   wortho = w;
      
#  Check orthogonality of qnew
#   test = norm(Q[:,1:k]'*qnew);
#   println("Orthogonality of new Vector: $test")

   return h,beta,wortho
end


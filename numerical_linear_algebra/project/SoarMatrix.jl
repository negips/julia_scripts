function SoarMatrix(ω,N)

#   using LinearAlgebra,ToeplitzMatrices

#  Inputs
#  N - No of Chebyshev points
#  ω - real angular frequency

#  Outputs - Matrices for the 2nd order Eigenvalue problem:
#  ( λ^2*[M] + λ*[D] + [K] )x = 0
#
# 
   M = 2;               # Second order derivates
   x,DM = chebdif(N,M);

#  Taking a case for suction boundary layer
   Re   =  300;     # Reynolds number (based on boundary layer thickness)
   δ    =  1.0;     # Boundary layer thickness
   U0   =  1.0;     # Free-stream velocity
   V0   = -1/Re;    # Wall suction velocity
   ν    =  1/Re;    # Kinematic velocity
   ymax =  1.;     # 100δ

#  eigenvector: [u v p]
# 
   Jac = (2. /ymax);
   DM[:,:,1] = DM[:,:,1].*Jac;
   DM[:,:,2] = DM[:,:,2].*(Jac^2);

#  Velocity profile
#  x goes from [1,-1]
#  Therefore the wall value is the last index not the first
   y = (x .- minimum(x)).*(ymax/2);
   U = (1. .- exp.(-y/δ));

#  Generate K
   eye = Matrix{Float64}(I,N,N);     # Identity
   zro = 0. .*eye;                   # Zero Matrix
   
   k1 = [(-ω*eye -(1/Re)*DM[:,:,2]) diagm(0 => DM[:,:,1]*U) zro];
   k2 = [zro  (-ω*eye -(1/Re)*DM[:,:,2]) DM[:,:,1]];
   k3 = [zro  DM[:,:,1]  zro];

   K  = vcat(k1, k2, k3);

   d1 = [diagm(0 => U) zro eye];
   d2 = [zro  diagm(0=> U) zro];
   d3 = [eye  zro  zro];

   D  = vcat(d1, d2, d3);

   m1 = [(-1/Re)*eye zro zro];
   m2 = [zro (-1/Re)*eye zro];
   m3 = [zro zro zro];
  
   M  = vcat(m1, m2, m3);

#  Remove Boundary conditions
#  Dirichlet at y=0,ymax
   i1  = collect(2:N-1);            # Dirichlet u_w=0., u_∞ =0.
   un  = length(i1); 

   i2  = collect(N+1:2*N-1);         # Dirichlet v_w =0., v_∞ =0.
   vn  = length(i2);

   i3  = collect(2*N+2:3*N);       # p∞ =0.    # arbitrary shift for pressure
   pn  = length(i3);

# 
   rws = vcat(i1,i2,i3);
   cls = vcat(i1,i2,i3);

   l1 = length(cls);
   z0 = Vector{Float64}(undef,l1)';
   z0 = 0. *z0;
#   println(cls)
#
   M = M[rws,cls];
   D = D[rws,cls];
   K = K[rws,cls];

#   M[end,:] = z0;
#   D[end,:] = z0;
#   K[end,:] = z0;

   BCk = hcat(zro,(-1/Re)*DM[:,:,2], DM[:,:,1] );
   k31 = BCk[1,cls];
   k3n = BCk[end,cls];

   K[un+vn+1,:] = k31;
#   K[end,:] = k3n;

#   K[end,un+vn+1:end] = DM[end,:,1];
#   K[end,un+1:un+vn]  = (-1/Re)*DM[end,i2 .- N,2];

   return x,y,U,DM,K,D,M,k1,k2,k3
end




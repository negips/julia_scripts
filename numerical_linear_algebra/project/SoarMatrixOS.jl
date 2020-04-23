function SoarMatrixOS(ω,β,N)

#   using LinearAlgebra,ToeplitzMatrices

#  Inputs
#  N - No of Chebyshev points
#  ω - real angular frequency

#  Outputs - Matrices for the 2nd order Eigenvalue problem:
#  ( λ^2*[M] + λ*[D] + [K] )x = 0
#
# 
   M = 4;               # Second order derivates
   x,DM = chebdif(N,M);

#  Taking a case for suction boundary layer
   Re   =  300;     # Reynolds number (based on boundary layer thickness)
   δ    =  1.0;     # Boundary layer thickness
   U0   =  1.0;     # Free-stream velocity
   V0   = -1/Re;    # Wall suction velocity
   ν    =  1/Re;    # Kinematic velocity
   ymax =  2.;     # 100δ

#  eigenvector: [u v p]
# 
   Jac = (2. /ymax);
   D1 = DM[:,:,1].*Jac;
   D2 = DM[:,:,2].*(Jac^2);
   D3 = DM[:,:,3].*(Jac^3);   
   D4 = DM[:,:,4].*(Jac^4);


#  Velocity profile
#  x goes from [1,-1]
#  Therefore the wall value is the last index not the first
   y = (x .- minimum(x)).*(ymax/2);
   U = (1. .- exp.(-y/δ));

#  Poiseuille flow
   Re   =  2000;     # Reynolds number (based on boundary layer thickness)
   y = (x).*(ymax/2);
   U = (1. .- y.^2);


#  Generate K
   eye = Matrix{Float64}(I,N,N);     # Identity
   zro = 0. .*eye;                   # Zero Matrix
   i = sqrt(Complex(-1.));                    # iota

#   tt = i*ω*(D2 - (β^2)*eye);
#   println(tt)
   
   k1 = [(i*ω*(D2 - (β^2)*eye) + (1/Re)*(D4 +(β^4)*eye -2*(β^2)*D2)) zro];
   k2 = [-i*β*(diagm(0=> D1*U)) (i*ω*eye + (1/Re)*(D2 - (β^2)*eye ))];

   K  = vcat(k1, k2);

   d1 = [(-i*diagm(0=> U)*(D2 - (β^2)*eye) + i*diagm(0=> D2*U) + (4/Re)*(β^2)*D1 - (4/Re)*D3) zro];
   d2 = [zro  (-i*diagm(0=> U) -(2/Re)*D1) ];

   D  = vcat(d1, d2);

   m1 = [(2*i*diagm(0=> U)*D1 + (4/Re)*D2) zro];
   m2 = [zro zro];
  
   M  = vcat(m1, m2);

#  Remove Boundary conditions
#  Dirichlet at y=0,ymax
   i1  = collect(2:N-1);            # Dirichlet u_w=0., u_∞ =0.
   un  = length(i1); 

   i2  = collect(N+2:2*N-1);        # Dirichlet v_w =0., v_∞ =0.
   vn  = length(i2);

   rws = vcat(i1,i2);
   cls = vcat(i1,i2);

#  Dv = 0. at the boundaries
#  => (D - αI)V = 0.

   l1 = length(cls);
   z0 = Vector{Float64}(undef,l1)';
   z0 = 0. *z0;

   M = M[rws,cls];
   D = D[rws,cls];
   K = K[rws,cls];
#
   M[1,:] = z0;
   D[1,:] = z0;
   K[1,:] = z0;

   M[un,:] = z0;
   D[un,:] = z0;
   K[un,:] = z0;

   BCk = hcat(D1, zro );
   BCd = hcat(-eye, zro );  
   k11 = BCk[1,cls];
   d11 = BCd[1,cls];

   k1n = BCk[end,cls];
   d1n = BCd[end,cls];
  
#   k1n = BCk[end,cls];

   K[1,:] = k11;
   D[1,:] = d11;

   K[un,:] = k1n;
   D[un,:] = d1n;

   return x,y,U,DM,K,D,M,k1,k2
end




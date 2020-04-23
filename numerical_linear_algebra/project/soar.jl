
#using Pkg
#Pkg.add("MATLAB")
using MatrixDepot,Random,TimerOutputs,PyPlot
using LinearAlgebra # for: eig, norm, etc
using Polynomials,MATLAB

const to = TimerOutput();

#include("arnoldi.jl")
include("GS_npass.jl")
#include("GS_modified.jl")

println("Second-Order Arnoldi (SOAR)")

# for Eigenvalue problems of the type:
# (λ^2M + λD + K)x = 0
# Kx = -λ(λM + D)

ngs = 1;    # No. of Gram-Schmidt passes

n=50;       # Matrix Size
kmax = 40;   # Max. Krylov dimension
if kmax>n
  kmax=n;
end  


# Diagonal M
M = Diagonal(fill(0.1, (n,n)));
Minv = Diagonal(fill(10., (n,n)));
D = Diagonal(fill(1.0, (n,n)));
d = fill(0.2,n);
ud = fill(-0.1,n-1);
ld = fill(-0.1,n-1);
K = Tridiagonal(ld, d, ud);
K[n,n] = 0.1;                 # modified boundary

# b = randn(Float64,n);
b = ones(n,1);


A = -Minv*D;
B = -Minv*K;

Q=fill(0.,(n,kmax));
Q2=fill(0.,(n,kmax));
P=fill(0.,(n,kmax));
T=fill(0.,(kmax+1,kmax) );

q = b/norm(b);
p = 0 .*q;
Q[:,1] = q;

#println(norm(b))

kk = 1;
Q2[:,kk] = q;           # non-zero vectors

for j in 1:kmax-1
  global q
  global p
  global r_ortho
  global tij
  global kk

  r1 = A*q;
  r2 = B*p;
  r = r1 + r2;
  s = q;
  tij,β,r_ortho = GS_npass(Q,r,j,ngs);      # Orthogonalize r with Q
  s_ortho = s - P[:,1:j]*tij;
  T[1:j,j]=tij;
  T[j+1,j]=β;

# Deflation  
  if abs(β)<10000*eps() || isnan(β)
    T[j+1,j] = 1.;
    β = 1.;
    q = fill(0., n);
    p = s_ortho;
  else    
    q = r_ortho/β;
    p = s_ortho/β;
  end 

  P[:,j+1]=p;
  Q[:,j+1]=q;
     
  if (norm(q)>10000*eps())
    kk=kk+1;
    Q2[:,kk]=q;
  end  

end   

kk=kk-1;
Mk = Q2[:,1:kk]'*M*Q2[:,1:kk];
Dk = Q2[:,1:kk]'*D*Q2[:,1:kk];
Kk = Q2[:,1:kk]'*K*Q2[:,1:kk];

m_Mk = mxarray(Mk);
m_Dk = mxarray(Dk);
m_Kk = mxarray(Kk);

mat"[$X,$E,$S]=polyeig($m_Kk,$m_Dk,$m_Mk)";

# d = diag(E);
e = sort(E);
plot(e, linestyle="none",marker=".");

# F = eigen(T[1:kmax,1:kmax]);
# 
# λr = real(F.values);
# λi = imag(F.values);
# 
# plot(λr,λi, linestyle="none",marker=".");

#to





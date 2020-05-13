println("Intializing SEM")

using PolynomialBases
using PyPlot,PyCall

close("all")

# define nodal bases
N           = 5 ;                               # polynomial degree
lx1         = N+1;                              # No of points
Basis       = LobattoLegendre(N)                # Polynomial Basis

Nd          = 8 ;                               # polynomial degree
lx1d        = Nd+1;                             # No of points
Basisd      = LobattoLegendre(Nd)               # Polynomial Basis

#basis2 = GaussLegendre(N)

xs          = 0.;                               # Domain start
xe          = 10.;                              # Domain end
nel         = 2;                                # No of elements
nnodes      = nel+1;                            # No of nodes
xc          = range(xs,stop=xe,length=nnodes);  # Element coordinates

zgm1        = Basis.nodes;                      # Reference GLL Points
wzm1        = Basis.weights;                    # Reference integration weights

zgm1d       = Basisd.nodes;                     # Reference GLL Points
wzm1d       = Basisd.weights;                   # Reference integration weights


# Reference Matrices
dxm1        = Basis.D;                          # Derivative Matrix
dxtm1       = Basis.D';                         # Derivative transpose

# GLL Points
xm1         = zeros(Float64,lx1,nel);
ym1         = zeros(Float64,lx1,nel);           # Tagging this along to make sense of derivatives
for i in 1:nel
  global xm1, ym1
 
  xi              = zeros(Float64,lx1);
  x0              = xc[i];
  x1              = xc[i+1];
  map_from_canonical!(xi,zgm1,x0,x1,Basis);
  xm1[:,i]        = xi;
  ym1[:,i]        = zgm1;
end 




# Local Geometric Matrices
xrm1   = zeros(Float64,lx1,nel);                # dx/dr
xsm1   = zeros(Float64,lx1,nel);                # dx/ds

yrm1   = zeros(Float64,lx1,nel);                # dy/dr
ysm1   = zeros(Float64,lx1,nel);                # dy/ds

rxm1   = zeros(Float64,lx1,nel);                # dr/dx
rym1   = zeros(Float64,lx1,nel);                # dr/dy
sxm1   = zeros(Float64,lx1,nel);                # ds/dx
sym1   = zeros(Float64,lx1,nel);                # ds/dy

jacm1  = zeros(Float64,lx1,nel);                # dr/dx
jacmi  = zeros(Float64,lx1,nel);                # dr/dx

for i in 1:nel
  global xrm1,xsm1,yrm1,ysm1 
  
  xrm1[:,i]  = dxm1*xm1[:,i];
  xsm1[:,i]  = 0 .*xm1[:,i];                    # Should be transpose

  yrm1[:,i]  = 0 .*ym1[:,i];
  ysm1[:,i]  = dxm1*ym1[:,i];                   # Should be transpose 

end


jacm1 = xrm1.*ysm1 - xsm1.*yrm1;
jacmi = 1 ./jacm1;                              # Inverse Jacobian

rxm1  = jacmi.*(ysm1);
sym1  = jacmi.*(xrm1);

rym1  = -jacmi.*(yrm1);
sxm1  = -jacmi.*(xsm1);

# Mass matrix
bm1   = jacm1.*wzm1;


# Gradient operator
gradx  = zeros(Float64,lx1,lx1,nel);            # d/dx
for i in 1:nel
  global gradx

  for j in 1:lx1
    gradx[j,:,i] = rxm1[j].*dxm1[j,:];
  end  
end

# Interpolation operator to de-aliased grid
intpm1d = zeros(Float64,lx1d,lx1);
intpm1d = interpolation_matrix(Basisd.nodes,Basis.nodes,Basis.baryweights);

jacm1d  = intpm1d*jacm1;
bm1d    = jacm1d.*wzm1d;            # Mass matrix on the de-aliased grid

gradxd = zeros(Float64,lx1d,lx1,nel);
bmd_matrix = zeros(Float64,lx1,lx1d,nel);

for i in 1:nel
  global gradxd    
  gradxd[:,:,i] = intpm1d*gradx[:,:,i];            # Interpolation of gradient on to the de-aliased grid
  vd            = intpm1d;
  bmd_matrix[:,:,i] = ones(Float64,lx1,1)*bm1d[:,i]';
#  bmd_matrix[:,:,i] = bmd_matrix[:,:,i]*vd;
end  

# Weak Laplacian
# wlp = -[ (BM1*∇v)^T.(∇) ]

wlp   = zeros(Float64,lx1,lx1,nel);
dvdx  = zeros(Float64,lx1,lx1,nel);

for i in 1:nel
  for j in 1:lx1
    global dvdx, wlp
    v    = zeros(Float64,lx1);
    v[j] = 1.;
    dv1  = gradx[:,:,i]*v;
    dvdx[:,j,i] = copy(dv1);
    dv1  = bm1[:,i].*dv1;
  
#    println(size(dv1))
#    println(dv1)
    wlp[j,:,i]  = -dv1'*gradx[:,:,i];
  end
end  


u = sin.(xm1[:,1]);
plot(xm1[:,1],u)

du = gradx[:,:,1]*u;
plot(xm1[:,1],du)


xd = intpm1d*xm1[:,1];
ud = intpm1d*u;
plot(xd,ud)

# # the function that will be interpolated
# ufunc(x) = sinpi(x); uprim(x) = π*cospi(x)
# #ufunc(x) = 1 / (1 + 25x^2); uprim(x) = -ufunc(x)^2*50x
# 
# for basis in (basis1, basis2)
#     u = ufunc.(basis.nodes)
# 
#     xplot = range(-1, stop=1, length=500)
#     uplot = interpolate(xplot, u, basis)
# 
#     fig1 = plot(xplot, ufunc.(xplot), label="u", xguide="x", yguide="u")
#     plot!(fig1, xplot, uplot, label="I(u)")
# 
#     fig2 = plot(xplot, uprim.(xplot), label="u'", xguide="x", yguide="u'")
#     plot!(fig2, xplot, interpolate(xplot, basis.D*u, basis), label="I(u)'")
# 
#     display(basis)
#     display(plot(fig1, fig2))
# end


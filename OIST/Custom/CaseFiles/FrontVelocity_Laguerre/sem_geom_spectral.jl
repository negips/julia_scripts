function sem_geom_laguerre(Basis,M2N,N2M,prec)
#     Generate the geometric matrices

#     Input:
#     N           : Degree of polynomial
#     Basis       : Structure for Polynomial basis of degree N
#     M2N         : Modal to Nodal Transform
#     N2M         : Nodal to Modal Transform

      VT    = prec
      nel   = 1

      lx1   = length(Basis.nodes);
      N     = lx1-1

      xm1   = zeros(VT,lx1);
      ym1   = zeros(VT,lx1);           # Tagging this along to make sense of derivatives
      xm1   = copy(Basis.nodes);
      ym1   = copy(Basis.nodes);
      wzm   = copy(Basis.weights);
       
      Wt    = exp.(-xm1)
     
#     Local Geometric Matrices
      xrm1   = zeros(VT,lx1);                # dx/dr
      xsm1   = zeros(VT,lx1);                # dx/ds
      
      yrm1   = zeros(VT,lx1);                # dy/dr
      ysm1   = zeros(VT,lx1);                # dy/ds
      
      rxm1   = zeros(VT,lx1);                # dr/dx
      rym1   = zeros(VT,lx1);                # dr/dy

      sxm1   = zeros(VT,lx1);                # ds/dx
      sym1   = zeros(VT,lx1);                # ds/dy
      
      jacm1  = zeros(VT,lx1);                # dr/dx
      jacmi  = zeros(VT,lx1);                # dr/dx
      
      dxm1   = Basis.D
      dxtm1  = Basis.D'
      for j in 1:lx1
        xrm1[j]  = 1.0
        xsm1[j]  = 0.0
      
        yrm1[j]  = 0.0
        ysm1[j]  = 1.0
      end  

      jacm1 = xrm1.*ysm1 - xsm1.*yrm1;
      jacmi = 1 ./jacm1;                              # Inverse Jacobian
      
      rxm1  = jacmi.*(ysm1);
      sym1  = jacmi.*(xrm1);
      
      rym1  = -jacmi.*(xsm1);
      sxm1  = -jacmi.*(yrm1);
      
      # Diagonal Mass matrix (as a vector)
      bm1   = Basis.weights.*Basis.orthoweights;
      
      # Gradient operator: Spectral -> Physical
      gradx  = copy(Basis.D)

      # Convective matrix assuming uniform velocity
      cnv    = (M2N')*diagm(Wt)*diagm(bm1)*gradx;
      
      # Weak Laplacian
      # wlp = -[ (BM1*∇v)^T.(∇) ]
      wlp  = -(gradx')*diagm(bm1)*dxm1

      # Laplacian (without integration by parts)      
      lap   = diagm(bm1)*gradx*gradx

      SpecGeom = SpectralMatrices(xm1,xrm1,rxm1,jacm1,jacmi,bm1,gradx,cnv,wlp,lap);

      println("SEM Spectral Geom: Done")

      return SpecGeom

end
#----------------------------------------------------------------------


function map_canonical_to_semiinfinite!(x,ξ,x0,scale;ifleft=false);

  l   = length(ξ)
  T   = eltype(x0)
  one = T(1)

  for i in 1:l
    if (ifleft)
      x[i] = x0 - scale*(one - ξ[i])/(one + ξ[i])
    else
      x[i] = x0 + scale*(one + ξ[i])/(one - ξ[i])
    end  
  end  

  return nothing
end

#---------------------------------------------------------------------- 
function map_semiinfinite_to_canonical!(ξ,x,x0,scale;ifleft=false);

  l   = length(x)
  T   = eltype(x)
  one = T(1)

  for i in 1:l
    if (ifleft)
      if (isinf(x[i]))
        ξ[i] = -one
      else  
        ξ[i] = (x0 - x[i] + scale)/(x0 - x[i] - scale)
      end  
    else
      if (isinf(x[i]))
        ξ[i] = one
      else
        ξ[i] = (x[i] - x0 - scale)/(x[i] - x0 + scale)
      end
    end  
  end  

  return nothing
end

#---------------------------------------------------------------------- 

struct SpectralMatrices
      xm1;
      xrm1;
      rxm1;
      jacm1;
      jacmi;
      bm1;
      gradx;
      cnv;
      wlp;
      lap;
end





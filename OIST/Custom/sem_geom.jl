function sem_geom(Basis,Basisd,xc,N,Nd,nel,dxm1,dxtm1,prec)
#     Generate the geometric matrices

#     Input:
#     N     : Degree of polynomial
#     Basis : Structure for Polynomial basis of degree N
#     Nd    : Degree of Dealiased polynomial
#     Basisd: Structure for Polynomial basis of degree Nd
#     nel   : No of elements
#     xc    : elemental nodes

      VT    = prec
#      VT = ComplexF64
#      VT    = Complex{prec}

      lx1  = N+1;
      lx1d = Nd+1;
#     GLL Points
      xm1         = zeros(VT,lx1,nel);
      ym1         = zeros(VT,lx1,nel);           # Tagging this along to make sense of derivatives
      zgm         = Basis.nodes;
      wzm         = Basis.weights;
      for i in 1:nel
        xi              = zeros(VT,lx1);
        x0              = xc[i];
        x1              = xc[i+1];
        map_from_canonical!(xi,zgm,x0,x1,Basis);
        xm1[:,i]        = xi;
        ym1[:,i]        = zgm;
      end 
      
#     Local Geometric Matrices
      xrm1   = zeros(VT,lx1,nel);                # dx/dr
      xsm1   = zeros(VT,lx1,nel);                # dx/ds
      
      yrm1   = zeros(VT,lx1,nel);                # dy/dr
      ysm1   = zeros(VT,lx1,nel);                # dy/ds
      
      rxm1   = zeros(VT,lx1,nel);                # dr/dx
      rym1   = zeros(VT,lx1,nel);                # dr/dy

      sxm1   = zeros(VT,lx1,nel);                # ds/dx
      sym1   = zeros(VT,lx1,nel);                # ds/dy
      
      jacm1  = zeros(VT,lx1,nel);                # dr/dx
      jacmi  = zeros(VT,lx1,nel);                # dr/dx
      
      for i in 1:nel
#        global xrm1,xsm1,yrm1,ysm1 
        
        xrm1[:,i]  = dxm1*xm1[:,i];
        xsm1[:,i]  = 0 .*xm1[:,i];                    # Should be transpose
      
        yrm1[:,i]  = 0 .*ym1[:,i];
        ysm1[:,i]  = dxm1*ym1[:,i];                   # Should be transpose 
      
      end
      
      
      jacm1 = xrm1.*ysm1 - xsm1.*yrm1;
      jacmi = 1 ./jacm1;                              # Inverse Jacobian
      
      rxm1  = jacmi.*(ysm1);
      sym1  = jacmi.*(xrm1);
      
      rym1  = -jacmi.*(xsm1);
      sxm1  = -jacmi.*(yrm1);
      
#     Diagonal Mass matrix (as a vector)
      bm1   = jacm1.*wzm;
      
#     Gradient operator
      gradx  = zeros(VT,lx1,lx1,nel);            # d/dx
      for i in 1:nel
        gradx[:,:,i]  = Diagonal(rxm1[:,i])*dxm1
      end  
      
#     Interpolation operator to de-aliased grid
      intpm1d = zeros(VT,lx1d,lx1);
      intpm1d = interpolation_matrix(Basisd.nodes,Basis.nodes,Basis.baryweights);
      
#     Matrices for Convection operator
      jacm1d  = intpm1d*jacm1;
      bm1d    = jacm1d.*Basisd.weights;            # Mass matrix on the de-aliased grid
      
      gradxd = zeros(VT,lx1d,lx1,nel);
      bmd_matrix = zeros(VT,lx1d,lx1d,nel);

      for i in 1:nel
        for j in 1:lx1d    
          bmd_matrix[j,j,i] = bm1d[j];
        end
      end
      
      bintpd = zeros(VT,lx1,lx1d,nel);        # Matrix to perform integration on the dealiased grid   
      for i in 1:nel
#       global gradxd,bmd_matrix 
        gradxd[:,:,i] = intpm1d*gradx[:,:,i];            # Interpolation of gradient on to the dealiased grid
        bintpd[:,:,i]   = (bmd_matrix[:,:,i]*intpm1d)';

      end

#     Convective matrix assuming uniform velocity
      cnv = zeros(VT,lx1,lx1,nel);
      for i in 1:nel
        cnv[:,:,i] = bintpd[:,:,i]*gradxd[:,:,i];
      end  
      
#     Weak Laplacian
#     wlp = -[ (BM1*∇v)^T.(∇) ]
      
      wlp   = zeros(VT,lx1,lx1,nel);
      dvdx  = zeros(VT,lx1,lx1,nel);
      
      for i in 1:nel
        dv1 = dxtm1 .+ 0.
        for j in 1:lx1
           dv1[j,:]  = bm1[:,i].*dxtm1[j,:].*rxm1[:,i] 
        end  
        dvdx[:,:,i] = dv1;
#        dv1  = bm1[:,i].*dv1;
        
        wlp[:,:,i]  = -dv1*gradx[:,:,i];
      end

#     Laplacian (without integration by parts)      
      lap   = zeros(VT,lx1,lx1,nel);
      
      for i in 1:nel
        lap[:,:,i] = Diagonal(bm1[:,i])*gradx[:,:,i]*gradx[:,:,i]
      end

      Geom = GeomMatrices(xm1,xrm1,rxm1,jacm1,jacmi,bm1,gradx,intpm1d,gradxd,bm1d,bintpd,cnv,wlp,lap,dvdx);

      println("SEM Geom: Done")

      return Geom

#      return xm1, xrm1, rxm1, jacm1, jacmi, bm1, gradx, intpm1d, gradxd, bm1d, wlp, dvdx 

end
#----------------------------------------------------------------------
function sem_geom(Basis,Basisd,xc::AbstractVector,N::Int,Nd::Int,nel::Int,dxm1,dxtm1,linf::Bool,rinf::Bool,Scale,prec)
#     Generate the geometric matrices

#     Input:
#     N     : Degree of polynomial
#     Basis : Structure for Polynomial basis of degree N
#     Nd    : Degree of Dealiased polynomial
#     Basisd: Structure for Polynomial basis of degree Nd
#     nel   : No of elements
#     xc    : elemental nodes
#     linf  : Left semi-infinite
#     rinf  : right semi-infinite

      VT    = prec
      one   = prec(1)
      zro   = prec(0)

      lx1  = N+1;
      lx1d = Nd+1;
#     GLL Points
      xm1         = zeros(VT,lx1,nel);
      ym1         = zeros(VT,lx1,nel);           # Tagging this along to make sense of derivatives
      zgm         = Basis.nodes;
      wzm         = Basis.weights;
      for i in 1:nel
        xi              = zeros(VT,lx1);
        x0              = xc[i];
        x1              = xc[i+1];
        if i==1 && linf
          map_canonical_to_semiinfinite!(xi,zgm,x1,Scale;ifleft=true)
        elseif i==nel && rinf
          map_canonical_to_semiinfinite!(xi,zgm,x0,Scale;ifleft=false)
        else
          map_from_canonical!(xi,zgm,x0,x1,Basis);
        end  
        xm1[:,i]        = xi;
        ym1[:,i]        = zgm;
      end 
      
#     Local Geometric Matrices
      xrm1   = zeros(VT,lx1,nel);                # dx/dr
      xsm1   = zeros(VT,lx1,nel);                # dx/ds
      
      yrm1   = zeros(VT,lx1,nel);                # dy/dr
      ysm1   = zeros(VT,lx1,nel);                # dy/ds
      
      rxm1   = zeros(VT,lx1,nel);                # dr/dx
      rym1   = zeros(VT,lx1,nel);                # dr/dy

      sxm1   = zeros(VT,lx1,nel);                # ds/dx
      sym1   = zeros(VT,lx1,nel);                # ds/dy
      
      jacm1  = zeros(VT,lx1,nel);                # dr/dx
      jacmi  = zeros(VT,lx1,nel);                # dr/dx
      
      for i in 1:nel
        # global xrm1,xsm1,yrm1,ysm1 
        if i==1 && linf
          for j in 1:lx1
            ξ           = Basis.nodes[j]
            xrm1[j,i]   = -(2.0*Scale)/((one + ξ)^2 )
            xsm1[j,i]   = zro
          end
        elseif i==nel && rinf
          for j in 1:lx1
            ξ           = Basis.nodes[j]
            xrm1[j,i]   = (2.0*Scale)/((one - ξ)^2 )
            xsm1[j,i]   = zro
          end
        else
          xrm1[:,i]  = dxm1*xm1[:,i];
          xsm1[:,i]  = 0 .*xm1[:,i];                  # Should be transpose
        end

        yrm1[:,i]  = 0 .*ym1[:,i];
        ysm1[:,i]  = dxm1*ym1[:,i];                   # Should be transpose 
      end         # i in 1:nel

      jacm1 = xrm1.*ysm1 - xsm1.*yrm1;
      jacmi = 1 ./jacm1;                              # Inverse Jacobian
      
      rxm1  = (ysm1);
      sym1  = (xrm1);
      
      rym1  = -(xsm1);
      sxm1  = -(yrm1);
      
#     Diagonal Mass matrix (as a vector)
      bm1   = jacm1.*wzm;
      
#     Gradient operator
      gradx  = zeros(VT,lx1,lx1,nel);            # d/dx
      for i in 1:nel
        gradx[:,:,i]  = Diagonal(jacmi[:,i].*rxm1[:,i])*dxm1
      end  
      
#     Interpolation operator to de-aliased grid
      intpm1d    = zeros(VT,lx1d,lx1);
      intpm1d    = interpolation_matrix(Basisd.nodes,Basis.nodes,Basis.baryweights);
      
#     Matrices for Convection operator
      jacm1d     = intpm1d*jacm1;
      bm1d       = jacm1d.*Basisd.weights;            # Mass matrix on the de-aliased grid
      gradxd     = zeros(VT,lx1d,lx1,nel);
      bmd_matrix = zeros(VT,lx1d,lx1d,nel);

      for i in 1:nel
        for j in 1:lx1d    
          bmd_matrix[j,j,i] = bm1d[j];
        end
      end
      
      bintpd = zeros(VT,lx1,lx1d,nel);        # Matrix to perform integration on the dealiased grid   
      for i in 1:nel
#       global gradxd,bmd_matrix 
        gradxd[:,:,i] = intpm1d*gradx[:,:,i];            # Interpolation of gradient on to the dealiased grid
        bintpd[:,:,i] = (bmd_matrix[:,:,i]*intpm1d)';
      end

#     Convective matrix assuming uniform velocity
      cnv = zeros(VT,lx1,lx1,nel);
      for i in 1:nel
        cnv[:,:,i] = diagm(Basis.weights)*diagm(rxm1[:,i])*dxm1
      end  
      
#     Weak Laplacian
#     wlp = -[ (BM1*∇v)^T.(∇) ]
      
      wlp   = zeros(VT,lx1,lx1,nel);
      dvdx  = zeros(VT,lx1,lx1,nel);
      for i in 1:nel
        # Canceling out one pair of Jacmi and Jac
        wlp[:,:,i] = -dxtm1*diagm(rxm1[:,i])*diagm(jacmi[:,i])'*diagm(Basis.weights)*diagm(rxm1[:,i])*dxm1
      end

#     Laplacian (without integration by parts) 
      lap   = zeros(VT,lx1,lx1,nel);
      for i in 1:nel
        lap[:,:,i] = diagm(Basis.weights)*diagm(rxm1[:,i])*dxm1*diagm(jacmi[:,i].*rxm1[:,i])*dxm1
      end

      Geom = GeomMatrices(xm1,xrm1,rxm1,jacm1,jacmi,bm1,gradx,intpm1d,gradxd,bm1d,bintpd,cnv,wlp,lap,dvdx);

      println("SEM Geom: Done")

      return Geom

#      return xm1, xrm1, rxm1, jacm1, jacmi, bm1, gradx, intpm1d, gradxd, bm1d, wlp, dvdx 

end
#---------------------------------------------------------------------- 
function sem_geom_laguerre(Basis,M2N,N2M,prec)
#     Generate the geometric matrices

#     Input:
#     N     : Degree of polynomial
#     Basis : Structure for Polynomial basis of degree N
#     Nd    : Degree of Dealiased polynomial
#     Basisd: Structure for Polynomial basis of degree Nd
#     nel   : No of elements
#     xc    : elemental nodes

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
      
#     Diagonal Mass matrix (as a vector)
      bm1   = Basis.weights;
      
#     Gradient operator
      gradx  = copy(Basis.D)
      
#     Interpolation operator to de-aliased grid
      intpm1d = zeros(VT,lx1,lx1);
      
#     Matrices for Convection operator
      jacm1d      = copy(jacm1);
      bm1d        = copy(Basis.weights);
      
      gradxd      = copy(Basis.D);
      bmd_matrix  = diagm(bm1d)
      
      bintpd      = zeros(VT,lx1,lx1);        # Matrix to perform integration on the dealiased grid   

#     Convective matrix assuming uniform velocity
      cnv         = zeros(VT,lx1,lx1);
      bintpd      = diagm(Basis.weights)
      cnv         = (M2N')*diagm(Wt)*bintpd*gradx;
      
#     Weak Laplacian
#     wlp = -[ (BM1*∇v)^T.(∇) ]
      
      dvdx  = zeros(VT,lx1,lx1);
      wlp  = -(gradx')*diagm(Wt.*Basis.weights)*dxm1

#     Laplacian (without integration by parts)      
      lap   = diagm(Wt.*bm1)*gradx*gradx

      Geom = GeomMatrices(xm1,xrm1,rxm1,jacm1,jacmi,bm1,gradx,intpm1d,gradxd,bm1d,bintpd,cnv,wlp,lap,dvdx);

      println("SEM Geom: Done")

      return Geom

#      return xm1, xrm1, rxm1, jacm1, jacmi, bm1, gradx, intpm1d, gradxd, bm1d, wlp, dvdx 

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


struct GeomMatrices
      xm1;
      xrm1;
      rxm1;
      jacm1;
      jacmi;
      bm1;
      gradx;
      intpm1d;
      gradxd;
      bm1d;
      bintpd;
      cnv;
      wlp;
      lap;
      dvdx;
end





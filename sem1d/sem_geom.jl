function sem_geom(Basis,Basisd,xc,N,Nd,nel,dxm1,dxtm1)
#     Generate the geometric matrices

#     Input:
#     N     : Degree of polynomial
#     Basis : Structure for Polynomial basis of degree N
#     Nd    : Degree of Dealiased polynomial
#     Basisd: Structure for Polynomial basis of degree Nd
#     nel   : No of elements
#     xc    : elemental nodes

      lx1  = N+1;
      lx1d = Nd+1;
#     GLL Points
      xm1         = zeros(Float64,lx1,nel);
      ym1         = zeros(Float64,lx1,nel);           # Tagging this along to make sense of derivatives
      for i in 1:nel
#        global xm1, ym1
       
        xi              = zeros(Float64,lx1);
        x0              = xc[i];
        x1              = xc[i+1];
        map_from_canonical!(xi,zgm1,x0,x1,Basis);
        xm1[:,i]        = xi;
        ym1[:,i]        = zgm1;
      end 
      
#     Local Geometric Matrices
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
#        global xrm1,xsm1,yrm1,ysm1 
        
        xrm1[:,i]  = dxm1*xm1[:,i];
        xsm1[:,i]  = 0 .*xm1[:,i];                    # Should be transpose
      
        yrm1[:,i]  = 0 .*ym1[:,i];
        ysm1[:,i]  = dxm1*ym1[:,i];                   # Should be transpose 
      
      end
      
      
      jacm1 = xrm1.*ysm1 - xsm1.*yrm1;
      jacmi = 1 ./jacm1;                              # Inverse Jacobian
      
      rxm1   = jacmi./xrm1

#      rxm1  = jacmi.*(ysm1);
#      sym1  = jacmi.*(xrm1);
      
#      rym1  = -jacmi.*(yrm1);
#      sxm1  = -jacmi.*(xsm1);
      
#     Diagonal Mass matrix (as a vector)
      bm1   = jacm1.*wzm1;
      
      
#     Gradient operator
      gradx  = zeros(Float64,lx1,lx1,nel);            # d/dx
      for i in 1:nel
        for j in 1:lx1
          gradx[j,:,i] = rxm1[j].*dxm1[j,:];
        end  
      end
      
#     Interpolation operator to de-aliased grid
      intpm1d = zeros(Float64,lx1d,lx1);
      intpm1d = interpolation_matrix(Basisd.nodes,Basis.nodes,Basis.baryweights);
      
#     Matrices for Convection operator
      jacm1d  = intpm1d*jacm1;
      tmp     = jacm1d.*wzm1d;            # Mass matrix on the de-aliased grid
      
      gradxd = zeros(Float64,lx1d,lx1,nel);
      bmd_matrix = zeros(Float64,lx1d,lx1d,nel);
      for i in 1:nel
        for j in 1:lx1d    
          bmd_matrix[j,j,i] = tmp[j];
        end
      end
      
      bm1d = zeros(Float64,lx1,lx1d,nel);        # Matrix to perform integration on the dealiased grid   
      for i in 1:nel
#        global gradxd,bmd_matrix 

        gradxd[:,:,i] = intpm1d*gradx[:,:,i];            # Interpolation of gradient on to the dealiased grid
        bm1d[:,:,i]   = (bmd_matrix[:,:,i]*intpm1d)';

      end

#     Convective matrix assuming uniform velocity
#      cnv = bm1d*gradxd;
      cnv = zeros(Float64,lx1,lx1,nel);
      for i in 1:nel
        cnv[:,:,i] = bm1d[:,:,i]*gradxd[:,:,i];
      end  
      
#     Weak Laplacian
#     wlp = -[ (BM1*∇v)^T.(∇) ]
      
      wlp   = zeros(Float64,lx1,lx1,nel);
      dvdx  = zeros(Float64,lx1,lx1,nel);
      
      for i in 1:nel
        dv1 = dxtm1 .+ 0.
        for j in 1:lx1
#          global dvdx, wlp
#          v    = zeros(Float64,lx1);
#          v[j] = 1.;
#          dv1  = gradx[:,:,i]*v;
#          println(size(rxm1))
#          println(size(dxtm1))
           dv1[j,:]  = bm1[:,i].*dxtm1[j,:].*rxm1[:,i] 
        end  
        dvdx[:,:,i] = copy(dv1);
#        dv1  = bm1[:,i].*dv1;
        
        wlp[:,:,i]  = -dv1*gradx[:,:,i];
      end
      
      Geom = GeomMatrices(xm1,xrm1,rxm1,jacm1,jacmi,bm1,gradx,intpm1d,gradxd,bm1d,cnv,wlp,dvdx);

      return Geom

#      return xm1, xrm1, rxm1, jacm1, jacmi, bm1, gradx, intpm1d, gradxd, bm1d, wlp, dvdx 

end

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
      cnv;
      wlp;
      dvdx;
end




function sem_geom2d(Bx,By,xc,yc,nel,prec)
#     Generate the geometric matrices

#     Input:
#     Bx    : Structure for Polynomial basis of degree N
#     By    : Structure for Polynomial basis of degree Nd
#     nel   : No of elements
#     xc    : x elemental nodes
#     yc    : y elemental nodes

      VT    = prec

      Nx    = PolynomialBases.degree(Bx)
      Ny    = PolynomialBases.degree(By)
      lx   = Nx+1;
      ly   = Ny+1;
#     GLL Points
      x           = zeros(VT,lx,ly,nel);
      y           = zeros(VT,lx,ly,nel);
      zgx         = Bx.nodes;
      wzx         = Bx.weights;
      zgy         = By.nodes;
      wzy         = By.weights;

      for i in 1:nel
#       global x, y
        xi              = zeros(VT,lx);
        x0              = xc[i];
        x1              = xc[i+1];
        map_from_canonical!(xi,zgx,x0,x1,Bx);

        yi              = zeros(VT,ly);
        y0              = yc[i];
        y1              = yc[i+1];
        map_from_canonical!(yi,zgy,y0,y1,By);
       
        x[:,:,i]        = kron(xi,ones(prec,1,ly))
        y[:,:,i]        = kron(ones(prec,lx,1),yi')
      end 
      
#     Local Geometric Matrices
      xr   = zeros(VT,lx,ly,nel);                # dx/dr
      xs   = zeros(VT,lx,ly,nel);                # dx/ds
      
      yr   = zeros(VT,lx,ly,nel);                # dy/dr
      ys   = zeros(VT,lx,ly,nel);                # dy/ds
      
      rx   = zeros(VT,lx,ly,nel);                # dr/dx
      ry   = zeros(VT,lx,ly,nel);                # dr/dy
      sx   = zeros(VT,lx,ly,nel);                # ds/dx
      sy   = zeros(VT,lx,ly,nel);                # ds/dy
      
      jac  = zeros(VT,lx,ly,nel);                # dr/dx
      jaci = zeros(VT,lx,ly,nel);                # dr/dx
#      
      for i in 1:nel
#        global xr,xs,yr,ys 
        
        xr[:,:,i]  = Bx.D*x[:,:,i];
        xs[:,:,i]  = x[:,:,i]*By.D';                    # Should be transpose
      
        yr[:,:,i]  = Bx.D*y[:,:,i];
        ys[:,:,i]  = y[:,:,i]*By.D';                   # Should be transpose 
      
      end
      
      det   = xr.*ys - xs.*yr;
      jac   = det
      jaci  = 1 ./jac;                              # Inverse Jacobian
      
      rx  = jaci.*(ys);
      sy  = jaci.*(xr);
      
      ry  = -jaci.*(xs);
      sx  = -jaci.*(yr);
      
#     Diagonal Mass matrix (as a vector)
      b   = jac.*kron(wzx,wzy');
#      
##     Gradient operator
#      gradx  = zeros(VT,lx,lx,nel);            # d/dx
#      for i in 1:nel
#        gradx[:,:,i]  = Diagonal(rx[:,i])*dx
#      end  
#      for i in 1:nel
#        for j in 1:lx
#          gradx[j,:,i] = rx[j,i].*dx[j,:];
#        end  
#      end
#      
##     Interpolation operator to de-aliased grid
#      intpd = zeros(VT,lxd,lx);
#      intpd = interpolation_matrix(Basisd.nodes,Basis.nodes,Basis.baryweights);
#      
##     Matrices for Convection operator
#      jacd  = intpd*jac;
#      bd    = jacd.*Basisd.weights;            # Mass matrix on the de-aliased grid
#      
#      gradxd = zeros(VT,lxd,lx,nel);
#      bmd_matrix = zeros(VT,lxd,lxd,nel);
#
#      for i in 1:nel
#        for j in 1:lxd    
#          bmd_matrix[j,j,i] = bd[j];
#        end
#      end
#      
#      bintpd = zeros(VT,lx,lxd,nel);        # Matrix to perform integration on the dealiased grid   
#      for i in 1:nel
##       global gradxd,bmd_matrix 
#        gradxd[:,:,i] = intpd*gradx[:,:,i];            # Interpolation of gradient on to the dealiased grid
#        bintpd[:,:,i]   = (bmd_matrix[:,:,i]*intpd)';
#
#      end
#
##     Convective matrix assuming uniform velocity
#      cnv = zeros(VT,lx,lx,nel);
#      for i in 1:nel
#        cnv[:,:,i] = bintpd[:,:,i]*gradxd[:,:,i];
#      end  
#      
##     Weak Laplacian
##     wlp = -[ (BM1*∇v)^T.(∇) ]
#      
#      wlp   = zeros(VT,lx,lx,nel);
#      dvdx  = zeros(VT,lx,lx,nel);
#      
#      for i in 1:nel
#        dv1 = dxt .+ 0.
#        for j in 1:lx
#           dv1[j,:]  = b[:,i].*dxt[j,:].*rx[:,i] 
#        end  
#        dvdx[:,:,i] = dv1;
##        dv1  = b[:,i].*dv1;
#        
#        wlp[:,:,i]  = -dv1*gradx[:,:,i];
#      end
      
      Geom = GeomMatrices2d(x,y,xr,xs,yr,ys,jac,jaci,rx,ry,sx,sy,b);

      println("SEM Geom: Done")

      return Geom

#      return x, xr, rx, jac, jacmi, b, gradx, intpd, gradxd, bd, wlp, dvdx 

end

mutable struct GeomMatrices2d
      x;
      y;
      xr;
      xs;
      yr;
      ys;
      jac;
      jaci;
      rx;
      ry;
      sx;
      sy;
      b;
#      rx;
#      jac;
#      jacmi;
#      b;
#      gradx;
#      intpd;
#      gradxd;
#      bd;
#      bintpd;
#      cnv;
#      wlp;
#      dvdx;
end





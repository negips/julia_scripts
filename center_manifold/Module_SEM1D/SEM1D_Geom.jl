function SEMGeoMat(Basis::BT,Basisd::BT,Inp::SEMInput;Dtype=Float64) where {BT<:NodalBasis} 
      # Generate the geometric matrices

      # Input:
      # N     : Degree of polynomial
      # Basis : Structure for Polynomial basis of degree N
      # Nd    : Degree of Dealiased polynomial
      # Basisd: Structure for Polynomial basis of degree Nd
      # nel   : No of elements
      # xc    : elemental nodes

      lx1  = Inp.lx1
      lx1d = Inp.lxd
      nel  = Inp.nel
      Dr   = Basis.D
      DrT  = Dr'

      # GLL/GL Points
      xm1         = zeros(Dtype,lx1,nel)
      zgm         = Basis.nodes
      for i in 1:nel
        xi              = zeros(Dtype,lx1)
        x0              = Inp.xc[i]
        x1              = Inp.xc[i+1]
        map_from_canonical!(xi,zgm,x0,x1,Basis)
        xm1[:,i]        = xi
      end 
      
      # Local Geometric Matrices
      xrm1   = zeros(Dtype,lx1,nel)                # dx/dr
      rxm1   = zeros(Dtype,lx1,nel)                # dr/dx
      
      jacm1  = zeros(Dtype,lx1,nel)                # dr/dx
      jacmi  = zeros(Dtype,lx1,nel)                # dr/dx
      
      for i in 1:nel
        xrm1[:,i]  = Dr*xm1[:,i]
      end
     
      # Jacobian
      jacm1 = copy(xrm1)
      jacmi = 1 ./jacm1                              # Inverse Jacobian
     
      # dr/dx
      rxm1  = copy(jacmi)
      
      # Diagonal Mass matrix (as a vector)
      bm1   = jacm1.*Basis.weights
      
      # Gradient operator
      Gradx  = zeros(Dtype,lx1,lx1,nel)            # d/dx
      for i in 1:nel
        Gradx[:,:,i]  = Diagonal(rxm1[:,i])*Dr
      end  
      
      # Interpolation operator to de-aliased grid
      Intpm1d = interpolation_matrix(Basisd.nodes,Basis.nodes,Basis.baryweights)
      
      # Matrices for Convection operator
      jacm1d  = Intpm1d*jacm1
      bm1d    = jacm1d.*Basisd.weights      # Mass (vector) on the de-aliased grid
      
      Gradxd  = zeros(Dtype,lx1d,lx1,nel)       # Interpolation of gradient on to the dealiased grid
      Bintpd  = zeros(Dtype,lx1,lx1d,nel)       # Matrix to perform integration on the dealiased grid   
      for i in 1:nel
        Gradxd[:,:,i]   = Intpm1d*Gradx[:,:,i]   
        Bintpd[:,:,i]   = (Diagonal(bm1d[:,i])*Intpm1d)'
      end

      # Dealiased Convective matrix assuming uniform velocity
      Convd = zeros(Dtype,lx1,lx1,nel)
      for i in 1:nel
        Convd[:,:,i] = Bintpd[:,:,i]*Gradxd[:,:,i]
      end  
      
      # Weak Laplacian
      # WLap = -[ (BM1*∇v)^T.(∇) ]
      WLap   = zeros(Dtype,lx1,lx1,nel)
      for i in 1:nel
        dvdx         = DrT*Diagonal(rxm1[:,i])*Diagonal(bm1[:,i])
        WLap[:,:,i]  = -dvdx*Gradx[:,:,i]
      end

      # Laplacian (without integration by parts)      
      Lap   = zeros(Dtype,lx1,lx1,nel)
      for i in 1:nel
        Lap[:,:,i] = Diagonal(bm1[:,i])*Gradx[:,:,i]*Gradx[:,:,i]
      end

      Geom = SEMGeoMat{Dtype}(xm1,xrm1,rxm1,jacm1,jacmi,bm1,Gradx,Intpm1d,Gradxd,bm1d,Bintpd,Convd,WLap,Lap)

      println("SEM Geom: Done")

      return Geom

end





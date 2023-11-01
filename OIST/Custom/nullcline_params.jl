#!/bin/julia
println("Setting Null-Cline Parameters")


if nullcline_set == 1
      nxA         = 1
      nyA         = 3
      nA          = nxA + nyA + 1      # No of free parameters/Conditions to satisfy
      
      FAC         = 0.2

      # Incident points     
      mA          = 3       
      xA          = zeros(Float64,mA)
      yA          = zeros(Float64,mA)
      if mA > 0
        xA        = FAC*[0.0 -2.0 -5.0]
        yA        = FAC*[0.0  0.7  6.0]
      end  
      
      # X-Derivative Points
      mxA = 0
      xdxA = zeros(Float64,mxA)
      ydxA = zeros(Float64,mxA)
      if (mx>0)
        xdxA = xA[2]
        ydxA = yA[2] 
      end  
      
      # Y-Derivative Points
      xdyA = zeros(Float64,mxA)
      ydyA = zeros(Float64,mxA)
      myA = 1
      if my>0
        xdyA[1] = xA[2]
        ydyA[1] = yA[2]
      end  
      
      fc = NullClineFcn(x,y,Xx,Yx,Xy,Yy,nx,ny)
      fc0 = fc[1]
      fcx = fc[2:nx+1]
      fcy = fc[nx+2:nx+ny+1]
#      
#      f(x,y) = FXY(x,y,fc0,fcx,fcy)
#      
#      xi = -20.0
#      yr0 = -50.0
#      yr1 = 50.0
#      dτ  = 1.0e-3
#      nsteps = 120000
#      
#      xx,yy = NullClines(f,xi,yr0,yr1,nsteps,dτ)
      
#      h1  = figure(num=1)
#      ax1 = h1.subplots()
#      
#      ax1.plot(xx,yy)
#      ax1.set_xlabel(L"x", fontsize=lafs)
#      ax1.set_ylabel(L"y", fontsize=lafs)
#      
#      ax1.plot(x,y,linestyle=" ",marker="x")
#      ax1.plot(Xx,Yx,linestyle=" ",marker="x")
#      ax1.plot(Xy,Yy,linestyle=" ",marker="x")
      
      # Function for G(x,y)
      #-------------------------------------------------- 
      nxB               = 1
      nyB               = 2
      nB                = nxB + nyB + 1      # No of free parameters/Conditions to satisfy
      
      # Incident points
      mB                = 2
      xB                = zeros(Float64,mB)
      yB                = zeros(Float64,mB)
      if mb > 0
        xB              = FAC*[0.0  13.0]
        yB              = FAC*[0.0  4.0]
      end  
      
      # X-Derivative Points
      mxB = 0
      xdxB              = zeros(Float64,mxB)
      ydxB              = zeros(Float64,mxB)
      if (mxB>0)
        xdxB            = xB[2]
        ydxB            = yB[2] 
      end  
      
      # Y-Derivative Points
      myB = 1
      xdyB              = zeros(Float64,myB)
      ydyB              = zeros(Float64,myB)
      
      if myB>0
        xdyB[1]         = xB[2]
        ydyB[1]         = yB[2]
      end  
      
      gc = NullClineFcn(x,y,Xx,Yx,Xy,Yy,nx,ny)
      gc0 = gc[1]
      gcx = gc[2:nx+1]
      gcy = gc[nx+2:nx+ny+1]
      
#      g(x,y) = FXY(x,y,gc0,gcx,gcy)
#      
#      xi = -10.0
#      yr0 = -50.0
#      yr1 = 1.0
#      dτ  = 1.0e-3
#      nsteps = 100000
#      
#      xx,yy = NullClines(g,xi,yr0,yr1,nsteps,dτ)



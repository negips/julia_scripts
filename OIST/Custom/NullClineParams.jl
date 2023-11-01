#!/bin/julia
struct NullClineParams
    
      nxA::Int
      nyA::Int
      nA::Int

      xA::Vector{Float64}
      yA::Vector{Float64}

      xdxA::Vector{Float64}
      ydxA::Vector{Float64}

      xdyA::Vector{Float64}
      ydyA::Vector{Float64}

      nxB::Int
      nyB::Int
      nB::Int

      xB::Vector{Float64}
      yB::Vector{Float64}

      xdxB::Vector{Float64}
      ydxB::Vector{Float64}

      xdyB::Vector{Float64}
      ydyB::Vector{Float64}

      fc0::Float64
      fcx::Vector{Float64}
      fcy::Vector{Float64}

      gc0::Float64
      gcx::Vector{Float64}
      gcy::Vector{Float64}

end
#---------------------------------------------------------------------- 
function GetNullClineXY(nullcline_set)

    if nullcline_set == 1
#     Foward-Backward Propagating pulses
      println("Parameters set for Forward-Backward Progagating pulses")

      nxA         = 1
      nyA         = 3
      nA          = nxA + nyA + 1      # No of free parameters/Conditions to satisfy
      
      FAC         = 0.2

      # Incident points     
      mA          = 3       
      xA          = zeros(Float64,mA)
      yA          = zeros(Float64,mA)
      if mA > 0
        xA        = FAC*[0.0; -2.0; -5.0]
        yA        = FAC*[0.0;  0.7;  6.0]
      end  
      
      # X-Derivative Points
      mxA = 0
      xdxA = zeros(Float64,mxA)
      ydxA = zeros(Float64,mxA)
      if (mxA>0)
        xdxA = xA[2]
        ydxA = yA[2] 
      end  
      
      # Y-Derivative Points
      myA = 1
      xdyA = zeros(Float64,myA)
      ydyA = zeros(Float64,myA)
      if myA>0
        xdyA[1] = xA[2]
        ydyA[1] = yA[2]
      end  
      
      # Points for G(x,y)
      #-------------------------------------------------- 
      nxB               = 1
      nyB               = 2
      nB                = nxB + nyB + 1      # No of free parameters/Conditions to satisfy
      
      # Incident points
      mB                = 2
      xB                = zeros(Float64,mB)
      yB                = zeros(Float64,mB)
      if mB > 0
        xB              = FAC*[0.0;  13.0]
        yB              = FAC*[0.0;  4.0]
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

    elseif nullcline_set == 2  
#     Foward-Backward Propagating pulses
      println("Parameters set for Limit Cycle Oscillation")

      nxA         = 1
      nyA         = 3
      nA          = nxA + nyA + 1      # No of free parameters/Conditions to satisfy
      
      FAC         = 0.25

      # Incident points     
      mA          = 3       
      xA          = zeros(Float64,mA)
      yA          = zeros(Float64,mA)
      if mA > 0
        xA        = FAC*[0.0; -10.0;  10.0]
        yA        = FAC*[0.0;  -5.0;   5.0]
      end  
      
      # X-Derivative Points
      mxA = 0
      xdxA = zeros(Float64,mxA)
      ydxA = zeros(Float64,mxA)
      if (mxA>0)
        xdxA = xA[2]
        ydxA = yA[2] 
      end  
      
      # Y-Derivative Points
      myA = 1
      xdyA = zeros(Float64,myA)
      ydyA = zeros(Float64,myA)
      if myA>0
        xdyA[1] = xA[2]
        ydyA[1] = yA[2]
      end  
      
      # Points for G(x,y)
      #-------------------------------------------------- 
      nxB               = 1
      nyB               = 2
      nB                = nxB + nyB + 1      # No of free parameters/Conditions to satisfy
      
      # Incident points
      mB                = 2
      xB                = zeros(Float64,mB)
      yB                = zeros(Float64,mB)
      if mB > 0
        xB              = FAC*[0.0;  25.0]
        yB              = FAC*[0.0;  8.0]
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


    else
      println("nullcline_set==$nullcline_set not defined")

      fc0 = []
      fcx = []
      fcy = []

      gc0 = []
      gcx = []
      gcy = []
     
    end  #  if nullcline_set == ?


    return nxA,nyA,nA,xA,yA,xdxA,ydxA,xdyA,ydyA,nxB,nyB,nB,xB,yB,xdxB,ydxB,xdyB,ydyB 
end

#---------------------------------------------------------------------- 
function GetNullClineParams(nullcline_set)

    nxA,nyA,nA,xA,yA,xdxA,ydxA,xdyA,ydyA,nxB,nyB,nB,xB,yB,xdxB,ydxB,xdyB,ydyB = GetNullClineXY(nullcline_set)
    fc = NullClineFcn(xA,yA,xdxA,ydxA,xdyA,ydyA,nxA,nyA)
    fc0 = fc[1]
    fcx = fc[2:nxA+1]
    fcy = fc[nxA+2:nxA+nyA+1]

    gc = NullClineFcn(xB,yB,xdxB,ydxB,xdyB,ydyB,nxB,nyB)
    gc0 = gc[1]
    gcx = gc[2:nxB+1]
    gcy = gc[nxB+2:nxB+nyB+1]

    params = NullClineParams(nxA,nyA,nA,xA,yA,xdxA,ydxA,xdyA,ydyA,nxB,nyB,nB,xB,yB,xdxB,ydxB,xdyB,ydyB,fc0,fcx,fcy,gc0,gcx,gcy)  

#    if nullcline_set == 1
##     Foward-Backward Propagating pulses
#      println("Parameters set for Forward-Backward Progagating pulses")
#
#      nxA         = 1
#      nyA         = 3
#      nA          = nxA + nyA + 1      # No of free parameters/Conditions to satisfy
#      
#      FAC         = 0.2
#
#      # Incident points     
#      mA          = 3       
#      xA          = zeros(Float64,mA)
#      yA          = zeros(Float64,mA)
#      if mA > 0
#        xA        = FAC*[0.0; -2.0; -5.0]
#        yA        = FAC*[0.0;  0.7;  6.0]
#      end  
#      
#      # X-Derivative Points
#      mxA = 0
#      xdxA = zeros(Float64,mxA)
#      ydxA = zeros(Float64,mxA)
#      if (mxA>0)
#        xdxA = xA[2]
#        ydxA = yA[2] 
#      end  
#      
#      # Y-Derivative Points
#      myA = 1
#      xdyA = zeros(Float64,myA)
#      ydyA = zeros(Float64,myA)
#      if myA>0
#        xdyA[1] = xA[2]
#        ydyA[1] = yA[2]
#      end  
#      
#      fc = NullClineFcn(xA,yA,xdxA,ydxA,xdyA,ydyA,nxA,nyA)
#      fc0 = fc[1]
#      fcx = fc[2:nxA+1]
#      fcy = fc[nxA+2:nxA+nyA+1]
#      
#      # Function for G(x,y)
#      #-------------------------------------------------- 
#      nxB               = 1
#      nyB               = 2
#      nB                = nxB + nyB + 1      # No of free parameters/Conditions to satisfy
#      
#      # Incident points
#      mB                = 2
#      xB                = zeros(Float64,mB)
#      yB                = zeros(Float64,mB)
#      if mB > 0
#        xB              = FAC*[0.0;  13.0]
#        yB              = FAC*[0.0;  4.0]
#      end  
#      
#      # X-Derivative Points
#      mxB = 0
#      xdxB              = zeros(Float64,mxB)
#      ydxB              = zeros(Float64,mxB)
#      if (mxB>0)
#        xdxB            = xB[2]
#        ydxB            = yB[2] 
#      end  
#      
#      # Y-Derivative Points
#      myB = 1
#      xdyB              = zeros(Float64,myB)
#      ydyB              = zeros(Float64,myB)
#      
#      if myB>0
#        xdyB[1]         = xB[2]
#        ydyB[1]         = yB[2]
#      end  
#      
#      gc = NullClineFcn(xB,yB,xdxB,ydxB,xdyB,ydyB,nxB,nyB)
#      gc0 = gc[1]
#      gcx = gc[2:nxB+1]
#      gcy = gc[nxB+2:nxB+nyB+1]
#
#    elseif
##     Foward-Backward Propagating pulses
#      println("Parameters set for Forward-Backward Progagating pulses")
#
#      nxA         = 1
#      nyA         = 3
#      nA          = nxA + nyA + 1      # No of free parameters/Conditions to satisfy
#      
#      FAC         = 0.2
#
#      # Incident points     
#      mA          = 3       
#      xA          = zeros(Float64,mA)
#      yA          = zeros(Float64,mA)
#      if mA > 0
#        xA        = FAC*[0.0; -2.0; -5.0]
#        yA        = FAC*[0.0;  0.7;  6.0]
#      end  
#      
#      # X-Derivative Points
#      mxA = 0
#      xdxA = zeros(Float64,mxA)
#      ydxA = zeros(Float64,mxA)
#      if (mxA>0)
#        xdxA = xA[2]
#        ydxA = yA[2] 
#      end  
#      
#      # Y-Derivative Points
#      myA = 1
#      xdyA = zeros(Float64,myA)
#      ydyA = zeros(Float64,myA)
#      if myA>0
#        xdyA[1] = xA[2]
#        ydyA[1] = yA[2]
#      end  
#      
#      fc = NullClineFcn(xA,yA,xdxA,ydxA,xdyA,ydyA,nxA,nyA)
#      fc0 = fc[1]
#      fcx = fc[2:nxA+1]
#      fcy = fc[nxA+2:nxA+nyA+1]
#      
#      # Function for G(x,y)
#      #-------------------------------------------------- 
#      nxB               = 1
#      nyB               = 2
#      nB                = nxB + nyB + 1      # No of free parameters/Conditions to satisfy
#      
#      # Incident points
#      mB                = 2
#      xB                = zeros(Float64,mB)
#      yB                = zeros(Float64,mB)
#      if mB > 0
#        xB              = FAC*[0.0;  13.0]
#        yB              = FAC*[0.0;  4.0]
#      end  
#      
#      # X-Derivative Points
#      mxB = 0
#      xdxB              = zeros(Float64,mxB)
#      ydxB              = zeros(Float64,mxB)
#      if (mxB>0)
#        xdxB            = xB[2]
#        ydxB            = yB[2] 
#      end  
#      
#      # Y-Derivative Points
#      myB = 1
#      xdyB              = zeros(Float64,myB)
#      ydyB              = zeros(Float64,myB)
#      
#      if myB>0
#        xdyB[1]         = xB[2]
#        ydyB[1]         = yB[2]
#      end  
#      
#      gc = NullClineFcn(xB,yB,xdxB,ydxB,xdyB,ydyB,nxB,nyB)
#      gc0 = gc[1]
#      gcx = gc[2:nxB+1]
#      gcy = gc[nxB+2:nxB+nyB+1]
#     
#
#    else
#      println("nullcline_set==$nullcline_set not defined")
#
#      fc0 = []
#      fcx = []
#      fcy = []
#
#      gc0 = []
#      gcx = []
#      gcy = []
#     
#    end  #  if nullcline_set == 1

    params = NullClineParams(nxA,nyA,nA,xA,yA,xdxA,ydxA,xdyA,ydyA,nxB,nyB,nB,xB,yB,xdxB,ydxB,xdyB,ydyB,fc0,fcx,fcy,gc0,gcx,gcy)  
    
    return params 
end
#---------------------------------------------------------------------- 


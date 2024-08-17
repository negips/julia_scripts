#!/bin/julia
mutable struct NullClineParams
    
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
@doc """

   function GetNullClineXY(nullcline_set)
   
   See GetNullClineParams for a list

"""
function GetNullClineXY(nullcline_set)

    if nullcline_set == 1
#----------------------------------------       
#     Foward-Backward Propagating pulses
      println("Forward-Backward Progagating pulses (puffs)")

      nxA         = 1
      nyA         = 3
      nA          = nxA + nyA + 1      # No of free parameters/Conditions to satisfy
      
      FAC         = 1.0

      # Incident points     
      mA          = 3       
      xA          = zeros(Float64,mA)
      yA          = zeros(Float64,mA)
      if mA > 0
#        xA        = FAC*[0.0; -0.4;    -1.0]
#        yA        = FAC*[0.0;  0.14;    1.2]
        xA        = FAC*[0.0; -0.4;    2.0]
        yA        = FAC*[0.0;  0.5;    5.0]
      end  
      
      # X-Derivative Points
      mxA = 0
      xdxA = zeros(Float64,mxA)
      ydxA = zeros(Float64,mxA)
      if (mxA>0)
        xdxA[1] = xA[2]
        ydxA[1] = yA[2] 
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
      #-------------------- 
      nxB               = 1
      nyB               = 2
      nB                = nxB + nyB + 1      # No of free parameters/Conditions to satisfy
      
      # Incident points
      mB                = 2
      xB                = zeros(Float64,mB)
      yB                = zeros(Float64,mB)
      if mB > 0
#        xB              = FAC*[0.0;  2.6]
#        yB              = FAC*[0.0;  0.8]
        xB              = FAC*[0.0;  5.5]
        yB              = FAC*[0.0;  4.0]
      end  
      
      # X-Derivative Points
      mxB = 0
      xdxB              = zeros(Float64,mxB)
      ydxB              = zeros(Float64,mxB)
      if (mxB>0)
        xdxB[1]            = xB[2]
        ydxB[1]            = yB[2] 
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
#----------------------------------------       
      println("Slugs. No Chnage in non-linear instability threshold")

      nxA         = 1
      nyA         = 3
      nA          = nxA + nyA + 1      # No of free parameters/Conditions to satisfy
      
      FAC         = 1.0

      # Incident points     
      mA          = 3       
      xA          = zeros(Float64,mA)
      yA          = zeros(Float64,mA)
      if mA > 0
#        xA        = FAC*[0.0; -0.4;    -1.0]
#        yA        = FAC*[0.0;  0.14;    1.2]

        xA        = FAC*[0.0; -0.4;    2.0]
        yA        = FAC*[0.0;  0.5;    5.0]
      end  
      
      # X-Derivative Points
      mxA = 0
      xdxA = zeros(Float64,mxA)
      ydxA = zeros(Float64,mxA)
      if (mxA>0)
        xdxA[1] = xA[2]
        ydxA[1] = yA[2] 
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
      #-------------------- 
      nxB               = 1
      nyB               = 2
      nB                = nxB + nyB + 1      # No of free parameters/Conditions to satisfy
      
      # Incident points
      mB                = 2
      xB                = zeros(Float64,mB)
      yB                = zeros(Float64,mB)
      if mB > 0
#        xB              = FAC*[0.0;  1.8]
#        yB              = FAC*[0.0;  0.8]

        xB              = FAC*[0.0;  3.5]
        yB              = FAC*[0.0;  5.0]
      end  
      
      # X-Derivative Points
      mxB = 0
      xdxB              = zeros(Float64,mxB)
      ydxB              = zeros(Float64,mxB)
      if (mxB>0)
        xdxB[1]            = xB[2]
        ydxB[1]            = yB[2] 
      end  
      
      # Y-Derivative Points
      myB = 1
      xdyB              = zeros(Float64,myB)
      ydyB              = zeros(Float64,myB)
      
      if myB>0
        xdyB[1]         = xB[2]
        ydyB[1]         = yB[2]
      end  

    elseif nullcline_set == 3
#----------------------------------------       
      println("Pulses. Smaller instability threshold.")

      nxA         = 1
      nyA         = 3
      nA          = nxA + nyA + 1      # No of free parameters/Conditions to satisfy
      
      FAC         = 1.0

      # Incident points     
      mA          = 3       
      xA          = zeros(Float64,mA)
      yA          = zeros(Float64,mA)
      if mA > 0
        xAshift   = [0.0; -0.25;    0.00]
        yAshift   = [0.0;  0.1883;  0.049]

        xA        = FAC*[0.0; -0.4;    2.0] .- xAshift
        yA        = FAC*[0.0;  0.5;    5.0] .- yAshift
      end  
      
      # X-Derivative Points
      mxA = 0
      xdxA = zeros(Float64,mxA)
      ydxA = zeros(Float64,mxA)
      if (mxA>0)
        xdxA[1] = xA[2]
        ydxA[1] = yA[2] 
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
      #-------------------- 
      nxB               = 1
      nyB               = 2
      nB                = nxB + nyB + 1      # No of free parameters/Conditions to satisfy
      
      # Incident points
      mB                = 2
      xB                = zeros(Float64,mB)
      yB                = zeros(Float64,mB)
      if mB > 0
#        xB              = FAC*[0.0;  1.8]
#        yB              = FAC*[0.0;  0.8]

        xB              = FAC*[0.0;  5.8]
        yB              = FAC*[0.0;  3.8]
      end  
      
      # X-Derivative Points
      mxB = 0
      xdxB              = zeros(Float64,mxB)
      ydxB              = zeros(Float64,mxB)
      if (mxB>0)
        xdxB[1]            = xB[2]
        ydxB[1]            = yB[2] 
      end  
      
      # Y-Derivative Points
      myB = 1
      xdyB              = zeros(Float64,myB)
      ydyB              = zeros(Float64,myB)
      
      if myB>0
        xdyB[1]         = xB[2]
        ydyB[1]         = yB[2]
      end  

    elseif nullcline_set == 4
#----------------------------------------       
      println("Slugs. Smaller instability threshold.")

      nxA         = 1
      nyA         = 3
      nA          = nxA + nyA + 1      # No of free parameters/Conditions to satisfy
      
      FAC         = 1.0

      # Incident points     
      mA          = 3       
      xA          = zeros(Float64,mA)
      yA          = zeros(Float64,mA)
      if mA > 0
        xAshift   = [0.0; -0.25;    0.00]
        yAshift   = [0.0;  0.1883;  0.049]

        xA        = FAC*[0.0; -0.4;    2.0] .- xAshift
        yA        = FAC*[0.0;  0.5;    5.0] .- yAshift
      end  
      
      # X-Derivative Points
      mxA = 0
      xdxA = zeros(Float64,mxA)
      ydxA = zeros(Float64,mxA)
      if (mxA>0)
        xdxA[1] = xA[2]
        ydxA[1] = yA[2] 
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
      #-------------------- 
      nxB               = 1
      nyB               = 2
      nB                = nxB + nyB + 1      # No of free parameters/Conditions to satisfy
      
      # Incident points
      mB                = 2
      xB                = zeros(Float64,mB)
      yB                = zeros(Float64,mB)
      if mB > 0
#        xB              = FAC*[0.0;  1.8]
#        yB              = FAC*[0.0;  0.8]

        xB              = FAC*[0.0;  3.6]
        yB              = FAC*[0.0;  3.8]
      end  
      
      # X-Derivative Points
      mxB = 0
      xdxB              = zeros(Float64,mxB)
      ydxB              = zeros(Float64,mxB)
      if (mxB>0)
        xdxB[1]            = xB[2]
        ydxB[1]            = yB[2] 
      end  
      
      # Y-Derivative Points
      myB = 1
      xdyB              = zeros(Float64,myB)
      ydyB              = zeros(Float64,myB)
      
      if myB>0
        xdyB[1]         = xB[2]
        ydyB[1]         = yB[2]
      end  
    elseif nullcline_set == 5
#----------------------------------------       
      println("Extreme Slugs. Smaller instability threshold.")

      nxA         = 1
      nyA         = 3
      nA          = nxA + nyA + 1      # No of free parameters/Conditions to satisfy
      
      FAC         = 1.0

      # Incident points     
      mA          = 3       
      xA          = zeros(Float64,mA)
      yA          = zeros(Float64,mA)
      if mA > 0
        xAshift   = [0.0; -0.25;    0.00]
        yAshift   = [0.0;  0.1883;  0.049]

        xA        = FAC*[0.0; -0.4;    2.0] .- xAshift
        yA        = FAC*[0.0;  0.5;    5.0] .- yAshift
      end  
      
      # X-Derivative Points
      mxA = 0
      xdxA = zeros(Float64,mxA)
      ydxA = zeros(Float64,mxA)
      if (mxA>0)
        xdxA[1] = xA[2]
        ydxA[1] = yA[2] 
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
      #-------------------- 
      nxB               = 1
      nyB               = 2
      nB                = nxB + nyB + 1      # No of free parameters/Conditions to satisfy
      
      # Incident points
      mB                = 2
      xB                = zeros(Float64,mB)
      yB                = zeros(Float64,mB)
      if mB > 0
#        xB              = FAC*[0.0;  1.8]
#        yB              = FAC*[0.0;  0.8]

        xB              = FAC*[0.0;  3.0]
        yB              = FAC*[0.0;  3.8]
      end  
      
      # X-Derivative Points
      mxB = 0
      xdxB              = zeros(Float64,mxB)
      ydxB              = zeros(Float64,mxB)
      if (mxB>0)
        xdxB[1]            = xB[2]
        ydxB[1]            = yB[2] 
      end  
      
      # Y-Derivative Points
      myB = 1
      xdyB              = zeros(Float64,myB)
      ydyB              = zeros(Float64,myB)
      
      if myB>0
        xdyB[1]         = xB[2]
        ydyB[1]         = yB[2]
      end 

    elseif (nullcline_set == 6)
#----------------------------------------       
#     Foward-Backward Propagating pulses
      println("Forward-Backward Progagating pulses (puffs)")

      nxA         = 1
      nyA         = 3
      nA          = nxA + nyA + 1      # No of free parameters/Conditions to satisfy
      
      FAC         = 1.0

      # Incident points     
      mA          = 3       
      xA          = zeros(Float64,mA)
      yA          = zeros(Float64,mA)
      if mA > 0
        xA        = FAC*[0.0; -0.4;    2.0]
        yA        = FAC*[0.0;  0.5;    4.0]
      end  
      
      # X-Derivative Points
      mxA = 0
      xdxA = zeros(Float64,mxA)
      ydxA = zeros(Float64,mxA)
      if (mxA>0)
        xdxA[1] = xA[2]
        ydxA[1] = yA[2] 
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
      #-------------------- 
      nxB               = 1
      nyB               = 2
      nB                = nxB + nyB + 1      # No of free parameters/Conditions to satisfy
      
      # Incident points
      mB                = 2
      xB                = zeros(Float64,mB)
      yB                = zeros(Float64,mB)
      if mB > 0
        xB              = FAC*[0.0;  7.5]
        yB              = FAC*[0.0;  4.0]
      end  
      
      # X-Derivative Points
      mxB = 0
      xdxB              = zeros(Float64,mxB)
      ydxB              = zeros(Float64,mxB)
      if (mxB>0)
        xdxB[1]            = xB[2]
        ydxB[1]            = yB[2] 
      end  
      
      # Y-Derivative Points
      myB = 1
      xdyB              = zeros(Float64,myB)
      ydyB              = zeros(Float64,myB)
      
      if myB>0
        xdyB[1]         = xB[2]
        ydyB[1]         = yB[2]
      end
    elseif nullcline_set == 7
#----------------------------------------       
      println("Parameters set for Extreme Slugs")

      nxA         = 1
      nyA         = 3
      nA          = nxA + nyA + 1      # No of free parameters/Conditions to satisfy
      
      FAC         = 1.0

      # Incident points     
      mA          = 3       
      xA          = zeros(Float64,mA)
      yA          = zeros(Float64,mA)
      if mA > 0
#        xAshift   = [0.0; -0.35;   -0.00]
#        yAshift   = [0.0;  0.0885;  0.016]
#
#        xA        = FAC*[0.0; -0.4;  -1.0] .- xAshift
#        yA        = FAC*[0.0;  0.14;  1.2] .- yAshift


        xAshift   = [0.0; -0.25;    0.00]
        yAshift   = [0.0;  0.1883;  0.049]

        xA        = FAC*[0.0; -0.4;    2.0] .- xAshift
        yA        = FAC*[0.0;  0.5;    5.0] .- yAshift
      end  
      
      # X-Derivative Points
      mxA = 0
      xdxA = zeros(Float64,mxA)
      ydxA = zeros(Float64,mxA)
      if (mxA>0)
        xdxA[1] = xA[2]
        ydxA[1] = yA[2] 
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
      #-------------------- 
      nxB               = 1
      nyB               = 2
      nB                = nxB + nyB + 1      # No of free parameters/Conditions to satisfy
      
      # Incident points
      mB                = 2
      xB                = zeros(Float64,mB)
      yB                = zeros(Float64,mB)
      if mB > 0
#        xB              = FAC*[0.0;  1.8]
#        yB              = FAC*[0.0;  0.8]

        xB              = FAC*[0.0;  3.5]
        yB              = FAC*[0.0;  5.0]
       
      end  
      
      # X-Derivative Points
      mxB = 0
      xdxB              = zeros(Float64,mxB)
      ydxB              = zeros(Float64,mxB)
      if (mxB>0)
        xdxB[1]            = xB[2]
        ydxB[1]            = yB[2] 
      end  
      
      # Y-Derivative Points
      myB = 1
      xdyB              = zeros(Float64,myB)
      ydyB              = zeros(Float64,myB)
      
      if myB>0
        xdyB[1]         = xB[2]
        ydyB[1]         = yB[2]
      end  
    elseif nullcline_set == 8
#----------------------------------------       
      println("Weak Slugs.")

      nxA         = 1
      nyA         = 3
      nA          = nxA + nyA + 1      # No of free parameters/Conditions to satisfy
      
      FAC         = 1.0

      # Incident points     
      mA          = 3       
      xA          = zeros(Float64,mA)
      yA          = zeros(Float64,mA)
      if mA > 0
#        xA        = FAC*[0.0; -0.4;    -1.0]
#        yA        = FAC*[0.0;  0.14;    1.2]

        xA        = FAC*[0.0; -0.4;    2.0]
        yA        = FAC*[0.0;  0.5;    5.0]
      end  
      
      # X-Derivative Points
      mxA = 0
      xdxA = zeros(Float64,mxA)
      ydxA = zeros(Float64,mxA)
      if (mxA>0)
        xdxA[1] = xA[2]
        ydxA[1] = yA[2] 
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
      #-------------------- 
      nxB               = 1
      nyB               = 2
      nB                = nxB + nyB + 1      # No of free parameters/Conditions to satisfy
      
      # Incident points
      mB                = 2
      xB                = zeros(Float64,mB)
      yB                = zeros(Float64,mB)
      if mB > 0
#        xB              = FAC*[0.0;  1.8]
#        yB              = FAC*[0.0;  0.8]

        xB              = FAC*[0.0;  4.0]
        yB              = FAC*[0.0;  5.0]
      end  
      
      # X-Derivative Points
      mxB = 0
      xdxB              = zeros(Float64,mxB)
      ydxB              = zeros(Float64,mxB)
      if (mxB>0)
        xdxB[1]            = xB[2]
        ydxB[1]            = yB[2] 
      end  
      
      # Y-Derivative Points
      myB = 1
      xdyB              = zeros(Float64,myB)
      ydyB              = zeros(Float64,myB)
      
      if myB>0
        xdyB[1]         = xB[2]
        ydyB[1]         = yB[2]
      end  

    elseif nullcline_set == 9
#----------------------------------------       
      println("Slugs. Linear Inhibitor")

      nxA         = 1
      nyA         = 3
      nA          = nxA + nyA + 1      # No of free parameters/Conditions to satisfy
      
      FAC         = 1.0

      # Incident points     
      mA          = 3       
      xA          = zeros(Float64,mA)
      yA          = zeros(Float64,mA)
      if mA > 0
        xA        = FAC*[0.0; -0.4;    2.0]
        yA        = FAC*[0.0;  0.5;    5.0]
      end  
      
      # X-Derivative Points
      mxA = 0
      xdxA = zeros(Float64,mxA)
      ydxA = zeros(Float64,mxA)
      if (mxA>0)
        xdxA[1] = xA[2]
        ydxA[1] = yA[2] 
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
      #-------------------- 
      nxB               = 1
      nyB               = 1
      nB                = nxB + nyB + 1      # No of free parameters/Conditions to satisfy
      
      # Incident points
      mB                = 2
      xB                = zeros(Float64,mB)
      yB                = zeros(Float64,mB)
      if mB > 0
        xB              = FAC*[0.0;  4.5]
        yB              = FAC*[0.0;  5.0]
      end  
      
      # X-Derivative Points
      mxB = 0
      xdxB              = zeros(Float64,mxB)
      ydxB              = zeros(Float64,mxB)
      if (mxB>0)
        xdxB[1]         = xB[2]
        ydxB[1]         = yB[2] 
      end  
      
      # Y-Derivative Points
      myB = 0
      xdyB              = zeros(Float64,myB)
      ydyB              = zeros(Float64,myB)
      
      if myB>0
        xdyB[1]         = xB[2]
        ydyB[1]         = yB[2]
      end  


    elseif nullcline_set == 10
#----------------------------------------       
      println("LCO. Activation Dominated")

      nxA         = 1
      nyA         = 3
      nA          = nxA + nyA + 1      # No of free parameters/Conditions to satisfy
      
      FAC         = 1.0

      # Incident points     
      mA          = 3       
      xA          = zeros(Float64,mA)
      yA          = zeros(Float64,mA)
      if mA > 0
#        xA        = FAC*[0.0; -10.0;  10.0]
#        yA        = FAC*[0.0;  -5.0;   5.0]

        xA        = FAC*[0.0; -0.4;   2.0]
        yA        = FAC*[0.0; -0.5;   4.0]
       
      end  
      
      # X-Derivative Points
      mxA = 0
      xdxA = zeros(Float64,mxA)
      ydxA = zeros(Float64,mxA)
      if (mxA>0)
        xdxA[1] = xA[2]
        ydxA[1] = yA[2] 
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
      #-------------------- 
      nxB               = 1
      nyB               = 3
      nB                = nxB + nyB + 1      # No of free parameters/Conditions to satisfy
      
      # Incident points
      mB                = 3
      xB                = zeros(Float64,mB)
      yB                = zeros(Float64,mB)
      if mB > 0
#        xB              = FAC*[0.0;  25.0]
#        yB              = FAC*[0.0;  8.0]
        xB              = FAC*[0.0;  7.0; -1.5]
        yB              = FAC*[0.0;  2.8; -0.4]
      end  
      
      # X-Derivative Points
      mxB = 0
      xdxB              = zeros(Float64,mxB)
      ydxB              = zeros(Float64,mxB)
      if (mxB>0)
        xdxB[1]            = xB[2]
        ydxB[1]            = yB[2] 
      end  
      
      # Y-Derivative Points
      myB = 1
      xdyB              = zeros(Float64,myB)
      ydyB              = zeros(Float64,myB)
      
      if myB>0
        xdyB[1]         = xB[2]
        ydyB[1]         = yB[2]
      end  

    elseif nullcline_set == 11
#----------------------------------------       
      println("LCO. Long Deactivation periods")

      nxA         = 1
      nyA         = 3
      nA          = nxA + nyA + 1      # No of free parameters/Conditions to satisfy
      
      FAC         = 1.0

      # Incident points     
      mA          = 3       
      xA          = zeros(Float64,mA)
      yA          = zeros(Float64,mA)
      if mA > 0
        xA        = FAC*[0.0; -0.4;  2.0]
        yA        = FAC*[0.0; -0.5;  4.0]
      end  
      
      # X-Derivative Points
      mxA = 0
      xdxA = zeros(Float64,mxA)
      ydxA = zeros(Float64,mxA)
      if (mxA>0)
        xdxA[1] = xA[2]
        ydxA[1] = yA[2] 
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
      #-------------------- 
      nxB               = 1
      nyB               = 3
      nB                = nxB + nyB + 1      # No of free parameters/Conditions to satisfy
      
      # Incident points
      mB                = 3
      xB                = zeros(Float64,mB)
      yB                = zeros(Float64,mB)
      if mB > 0
#        xB              = FAC*[0.0; -0.55; 7.0]
#        yB              = FAC*[0.0; -0.60; 3.0]
        xB              = FAC*[0.0; -0.55; 9.0]
        yB              = FAC*[0.0; -0.60; 3.0]
      end  
      
      # X-Derivative Points
      mxB = 0
      xdxB              = zeros(Float64,mxB)
      ydxB              = zeros(Float64,mxB)
      if (mxB>0)
        xdxB[1]            = xB[2]
        ydxB[1]            = yB[2] 
      end  
      
      # Y-Derivative Points
      myB = 1
      xdyB              = zeros(Float64,myB)
      ydyB              = zeros(Float64,myB)
      
      if myB>0
        xdyB[1]         = xB[2]
        ydyB[1]         = yB[2]
      end  

    elseif nullcline_set == 12
#----------------------------------------       
      println("LCO. Long Deactivation periods")

      nxA         = 1
      nyA         = 3
      nA          = nxA + nyA + 1      # No of free parameters/Conditions to satisfy
      
      FAC         = 1.0

      # Incident points     
      mA          = 3       
      xA          = zeros(Float64,mA)
      yA          = zeros(Float64,mA)
      if mA > 0
        xA        = FAC*[0.0; -0.4;    2.0]
        yA        = FAC*[0.0; -0.5;    4.0]
      end  
      
      # X-Derivative Points
      mxA = 0
      xdxA = zeros(Float64,mxA)
      ydxA = zeros(Float64,mxA)
      if (mxA>0)
        xdxA[1] = xA[2]
        ydxA[1] = yA[2] 
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
      #-------------------- 
      nxB               = 3
      nyB               = 1
      nB                = nxB + nyB + 1      # No of free parameters/Conditions to satisfy
      
      # Incident points
      mB                = 4
      xB                = zeros(Float64,mB)
      yB                = zeros(Float64,mB)
      if mB > 0
        xB              = FAC*[0.0; -1.00; 6.5; -3.0]
        yB              = FAC*[0.0; -0.60; 3.0; -3.0]
      end  
      
      # X-Derivative Points
      mxB = 0
      xdxB              = zeros(Float64,mxB)
      ydxB              = zeros(Float64,mxB)
      if (mxB>0)
        xdxB[1]            = xB[1]
        ydxB[1]            = yB[1] 
      end  
      
      # Y-Derivative Points
      myB = 0
      xdyB              = zeros(Float64,myB)
      ydyB              = zeros(Float64,myB)
      
      if myB>0
        xdyB[1]         = xB[2]
        ydyB[1]         = yB[2]
      end  

    elseif nullcline_set == 13
#----------------------------------------       
      println("Two fixed Point on the upper and lower branches")

      nxA         = 1
      nyA         = 3
      nA          = nxA + nyA + 1      # No of free parameters/Conditions to satisfy
      
      FAC         = 1.0

      # Incident points     
      mA          = 3       
      xA          = zeros(Float64,mA)
      yA          = zeros(Float64,mA)
      if mA > 0
        xA        = FAC*[0.0;  -3.0;  2.0]
        yA        = FAC*[0.0;  -3.0;  4.0]
      end  
      
      # X-Derivative Points
      mxA = 0
      xdxA = zeros(Float64,mxA)
      ydxA = zeros(Float64,mxA)
      if (mxA>0)
        xdxA[1] = xA[2]
        ydxA[1] = yA[2] 
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
      #-------------------- 
      nxB               = 1
      nyB               = 1
      nB                = nxB + nyB + 1      # No of free parameters/Conditions to satisfy
      
      # Incident points
      mB                = 2
      xB                = zeros(Float64,mB)
      yB                = zeros(Float64,mB)
      if mB > 0
        xB              = FAC*[0.0;  3.0]
        yB              = FAC*[0.0; -6.0]
       
      end  
      
      # X-Derivative Points
      mxB = 0
      xdxB              = zeros(Float64,mxB)
      ydxB              = zeros(Float64,mxB)
      if (mxB>0)
        xdxB[1]            = xB[2]
        ydxB[1]            = yB[2] 
      end  
      
      # Y-Derivative Points
      myB = 0
      xdyB              = zeros(Float64,myB)
      ydyB              = zeros(Float64,myB)
      
      if myB>0
        xdyB[1]         = xB[2]
        ydyB[1]         = yB[2]
      end  

    elseif nullcline_set == 14
#----------------------------------------       
      println("Symmetric fixed points on the upper and lower branches")

      nxA         = 1
      nyA         = 3
      nA          = nxA + nyA + 1      # No of free parameters/Conditions to satisfy
      
      FAC         = 1.0

      # Incident points     
      mA          = 3       
      xA          = zeros(Float64,mA)
      yA          = zeros(Float64,mA)
      if mA > 0
        xA        = FAC*[0.0;  -3.0;  3.0]
        yA        = FAC*[0.0;  -3.0;  3.0]
      end  
      
      # X-Derivative Points
      mxA = 0
      xdxA = zeros(Float64,mxA)
      ydxA = zeros(Float64,mxA)
      if (mxA>0)
        xdxA[1] = xA[2]
        ydxA[1] = yA[2] 
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
      #-------------------- 
      nxB               = 1
      nyB               = 1
      nB                = nxB + nyB + 1      # No of free parameters/Conditions to satisfy
      
      # Incident points
      mB                = 2
      xB                = zeros(Float64,mB)
      yB                = zeros(Float64,mB)
      if mB > 0
        xB              = FAC*[0.0;  2.0]
        yB              = FAC*[0.0; -8.0]
       
      end  
      
      # X-Derivative Points
      mxB = 0
      xdxB              = zeros(Float64,mxB)
      ydxB              = zeros(Float64,mxB)
      if (mxB>0)
        xdxB[1]            = xB[2]
        ydxB[1]            = yB[2] 
      end  
      
      # Y-Derivative Points
      myB = 0
      xdyB              = zeros(Float64,myB)
      ydyB              = zeros(Float64,myB)
      
      if myB>0
        xdyB[1]         = xB[2]
        ydyB[1]         = yB[2]
      end

    elseif nullcline_set == 15
#----------------------------------------       
      println("Symmetric LCO")

      nxA         = 1
      nyA         = 3
      nA          = nxA + nyA + 1      # No of free parameters/Conditions to satisfy
      
      FAC         = 1.0

      # Incident points     
      mA          = 3       
      xA          = zeros(Float64,mA)
      yA          = zeros(Float64,mA)
      if mA > 0
        xA        = FAC*[0.0;  -3.0;  3.0]
        yA        = FAC*[0.0;  -3.0;  3.0]
      end  
      
      # X-Derivative Points
      mxA = 0
      xdxA = zeros(Float64,mxA)
      ydxA = zeros(Float64,mxA)
      if (mxA>0)
        xdxA[1] = xA[2]
        ydxA[1] = yA[2] 
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
      #-------------------- 
      nxB               = 1
      nyB               = 1
      nB                = nxB + nyB + 1      # No of free parameters/Conditions to satisfy
      
      # Incident points
      mB                = 2
      xB                = zeros(Float64,mB)
      yB                = zeros(Float64,mB)
      if mB > 0
        xB              = FAC*[0.0;  6.0]
        yB              = FAC*[0.0;  2.0]
       
      end  
      
      # X-Derivative Points
      mxB = 0
      xdxB              = zeros(Float64,mxB)
      ydxB              = zeros(Float64,mxB)
      if (mxB>0)
        xdxB[1]            = xB[2]
        ydxB[1]            = yB[2] 
      end  
      
      # Y-Derivative Points
      myB = 0
      xdyB              = zeros(Float64,myB)
      ydyB              = zeros(Float64,myB)
      
      if myB>0
        xdyB[1]         = xB[2]
        ydyB[1]         = yB[2]
      end  

    elseif nullcline_set == 16
#----------------------------------------       
      println("Asymmetric fixed points")

      nxA         = 1
      nyA         = 3
      nA          = nxA + nyA + 1      # No of free parameters/Conditions to satisfy
      
      FAC         = 1.0

      # Incident points     
      mA          = 3       
      xA          = zeros(Float64,mA)
      yA          = zeros(Float64,mA)
      if mA > 0
        xA        = FAC*[0.0; -0.4;    2.0]
        yA        = FAC*[0.0;  0.5;    5.0]
      end  
      
      # X-Derivative Points
      mxA = 0
      xdxA = zeros(Float64,mxA)
      ydxA = zeros(Float64,mxA)
      if (mxA>0)
        xdxA[1] = xA[2]
        ydxA[1] = yA[2] 
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
      #-------------------- 
      nxB               = 1
      nyB               = 2
      nB                = nxB + nyB + 1      # No of free parameters/Conditions to satisfy
      
      # Incident points
      mB                = 2
      xB                = zeros(Float64,mB)
      yB                = zeros(Float64,mB)
      if mB > 0
        xB              = FAC*[0.0;  3.5]
        yB              = FAC*[0.0;  3.0]
      end  
      
      # X-Derivative Points
      mxB = 0
      xdxB              = zeros(Float64,mxB)
      ydxB              = zeros(Float64,mxB)
      if (mxB>0)
        xdxB[1]            = xB[2]
        ydxB[1]            = yB[2] 
      end  
      
      # Y-Derivative Points
      myB = 1
      xdyB              = zeros(Float64,myB)
      ydyB              = zeros(Float64,myB)
      
      if myB>0
        xdyB[1]         = xB[2]
        ydyB[1]         = yB[2]
      end  
    elseif nullcline_set == 17
#----------------------------------------       
      println("Unstable G-nullcline")

      nxA         = 1
      nyA         = 3
      nA          = nxA + nyA + 1      # No of free parameters/Conditions to satisfy
      
      FAC         = 1.0

      # Incident points     
      mA          = 3       
      xA          = zeros(Float64,mA)
      yA          = zeros(Float64,mA)
      if mA > 0
        xA        = FAC*[0.0;  -3.0;  3.0]
        yA        = FAC*[0.0;  -3.0;  3.0]
       end  
      
      # X-Derivative Points
      mxA = 0
      xdxA = zeros(Float64,mxA)
      ydxA = zeros(Float64,mxA)
      if (mxA>0)
        xdxA[1] = xA[2]
        ydxA[1] = yA[2] 
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
      #-------------------- 
      nxB               = 1
      nyB               = 1
      nB                = nxB + nyB + 1      # No of free parameters/Conditions to satisfy
      
      # Incident points
      mB                = 2
      xB                = zeros(Float64,mB)
      yB                = zeros(Float64,mB)
      if mB > 0
        xB              = FAC*[0.0;  10.0]
        yB              = FAC*[0.0;  -5.0]
      end  
      
      # X-Derivative Points
      mxB = 0
      xdxB              = zeros(Float64,mxB)
      ydxB              = zeros(Float64,mxB)
      if (mxB>0)
        xdxB[1]            = xB[2]
        ydxB[1]            = yB[2] 
      end  
      
      # Y-Derivative Points
      myB = 0
      xdyB              = zeros(Float64,myB)
      ydyB              = zeros(Float64,myB)
      
      if myB>0
        xdyB[1]         = xB[2]
        ydyB[1]         = yB[2]
      end  

    elseif nullcline_set == 18
#----------------------------------------       
      println("Marginal Instability")

      nxA         = 1
      nyA         = 3
      nA          = nxA + nyA + 1      # No of free parameters/Conditions to satisfy
      
      FAC         = 1.0

      # Incident points     
      mA          = 3 
      xA          = zeros(Float64,mA)
      yA          = zeros(Float64,mA)
      if mA > 0
        xA        = FAC*[0.0;   0.0;   2.0]
        yA        = FAC*[0.0;   6.0;   4.5]
       end  
      
      # X-Derivative Points
      mxA = 0
      xdxA = zeros(Float64,mxA)
      ydxA = zeros(Float64,mxA)
      if (mxA>0)
        xdxA[1] = xA[2]
        ydxA[1] = yA[2] 
      end  
      
      # Y-Derivative Points
      myA = 1
      xdyA = zeros(Float64,myA)
      
      ydyA = zeros(Float64,myA)
      if myA>0
        xdyA[1] = xA[1]
        ydyA[1] = yA[1]
      end  
      
      # Points for G(x,y)
      #-------------------- 
      nxB               = 1
      nyB               = 1
      nB                = nxB + nyB + 1      # No of free parameters/Conditions to satisfy
      
      # Incident points
      mB                = 2
      xB                = zeros(Float64,mB)
      yB                = zeros(Float64,mB)
      if mB > 0
        xB              = FAC*[0.0;  2.6]
        yB              = FAC*[0.0;  4.0]
      end  
      
      # X-Derivative Points
      mxB = 0
      xdxB              = zeros(Float64,mxB)
      ydxB              = zeros(Float64,mxB)
      if (mxB>0)
        xdxB[1]            = xB[2]
        ydxB[1]            = yB[2] 
      end  
      
      # Y-Derivative Points
      myB = 0
      xdyB              = zeros(Float64,myB)
      ydyB              = zeros(Float64,myB)
      
      if myB>0
        xdyB[1]         = xB[2]
        ydyB[1]         = yB[2]
      end  

    elseif nullcline_set == 19
#----------------------------------------       
#     Stable oscillatory linear mode
      println("Stable Oscillatory Linear Mode")

      nxA         = 1
      nyA         = 3
      nA          = nxA + nyA + 1      # No of free parameters/Conditions to satisfy
      
      FAC         = 1.0

      # Incident points     
      mA          = 3       
      xA          = zeros(Float64,mA)
      yA          = zeros(Float64,mA)
      if mA > 0
#        xA        = FAC*[0.0; -0.4;    -1.0]
#        yA        = FAC*[0.0;  0.14;    1.2]
        xA        = FAC*[0.0; -0.4;    2.0]
        yA        = FAC*[0.0;  0.5;    5.0]
#        xA        = FAC*[0.0;   0.0;   -0.001]
#        yA        = FAC*[0.0;   6.0;   0.05]
      end  
      
      # X-Derivative Points
      mxA = 0
      xdxA = zeros(Float64,mxA)
      ydxA = zeros(Float64,mxA)
      if (mxA>0)
        xdxA[1] = xA[2]
        ydxA[1] = yA[2] 
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
      #-------------------- 
      nxB               = 2
      nyB               = 1
      nB                = nxB + nyB + 1      # No of free parameters/Conditions to satisfy
      
      # Incident points
      mB                = 2
      xB                = zeros(Float64,mB)
      yB                = zeros(Float64,mB)
      if mB > 0
#        xB              = FAC*[0.0;  2.6]
#        yB              = FAC*[0.0;  0.8]
        xB              = FAC*[0.0;  5.0]
        yB              = FAC*[0.0;  6.0]
      end  
      
      # X-Derivative Points
      mxB = 1
      xdxB              = zeros(Float64,mxB)
      ydxB              = zeros(Float64,mxB)
      if (mxB>0)
        xdxB[1]            = xB[2]
        ydxB[1]            = yB[2] 
      end  
      
      # Y-Derivative Points
      myB = 0
      xdyB              = zeros(Float64,myB)
      ydyB              = zeros(Float64,myB)
      
      if myB>0
        xdyB[1]         = xB[2]
        ydyB[1]         = yB[2]
      end

    elseif nullcline_set == 20
#----------------------------------------       
      println("Slugs. Large non-linear instability threshold")

      nxA         = 1
      nyA         = 3
      nA          = nxA + nyA + 1      # No of free parameters/Conditions to satisfy
      
      FAC         = 1.0

      # Incident points     
      mA          = 3       
      xA          = zeros(Float64,mA)
      yA          = zeros(Float64,mA)
      if mA > 0
        xA        = FAC*[0.0; -0.5;    2.0]
        yA        = FAC*[0.0;  0.6;    5.0]
      end  
      
      # X-Derivative Points
      mxA = 0
      xdxA = zeros(Float64,mxA)
      ydxA = zeros(Float64,mxA)
      if (mxA>0)
        xdxA[1] = xA[2]
        ydxA[1] = yA[2] 
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
      #-------------------- 
      nxB               = 1
      nyB               = 1
      nB                = nxB + nyB + 1      # No of free parameters/Conditions to satisfy
      
      # Incident points
      mB                = 2
      xB                = zeros(Float64,mB)
      yB                = zeros(Float64,mB)
      if mB > 0
        xB              = FAC*[0.0;  3.5]
        yB              = FAC*[0.0;  4.6]
      end 
      
      # X-Derivative Points
      mxB = 0
      xdxB              = zeros(Float64,mxB)
      ydxB              = zeros(Float64,mxB)
      if (mxB>0)
        xdxB[1]            = xB[2]
        ydxB[1]            = yB[2] 
      end  
      
      # Y-Derivative Points
      myB = 0
      xdyB              = zeros(Float64,myB)
      ydyB              = zeros(Float64,myB)
      
      if myB>0
        xdyB[1]         = xB[2]
        ydyB[1]         = yB[2]
      end  

    elseif nullcline_set == 21
#----------------------------------------       
      println("Slugs. Small non-linear instability threshold")

      nxA         = 1
      nyA         = 3
      nA          = nxA + nyA + 1      # No of free parameters/Conditions to satisfy
      
      FAC         = 1.0

      # Incident points     
      mA          = 3       
      xA          = zeros(Float64,mA)
      yA          = zeros(Float64,mA)
      if mA > 0
        xA        = FAC*[0.0; -0.15;   2.0]
        yA        = FAC*[0.0;  0.32;   4.0]
      end  
      
      # X-Derivative Points
      mxA = 0
      xdxA = zeros(Float64,mxA)
      ydxA = zeros(Float64,mxA)
      if (mxA>0)
        xdxA[1] = xA[2]
        ydxA[1] = yA[2] 
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
      #-------------------- 
      nxB               = 1
      nyB               = 1
      nB                = nxB + nyB + 1      # No of free parameters/Conditions to satisfy
      
      # Incident points
      mB                = 2
      xB                = zeros(Float64,mB)
      yB                = zeros(Float64,mB)
      if mB > 0
        xB              = FAC*[0.0;  4.0]
        yB              = FAC*[0.0;  4.0]
      end 
      
      # X-Derivative Points
      mxB = 0
      xdxB              = zeros(Float64,mxB)
      ydxB              = zeros(Float64,mxB)
      if (mxB>0)
        xdxB[1]            = xB[2]
        ydxB[1]            = yB[2] 
      end  
      
      # Y-Derivative Points
      myB = 0
      xdyB              = zeros(Float64,myB)
      ydyB              = zeros(Float64,myB)
      
      if myB>0
        xdyB[1]         = xB[2]
        ydyB[1]         = yB[2]
      end  

    elseif nullcline_set == 22
#----------------------------------------       
      println("Slugs. Small non-linear instability threshold (2)")

      nxA         = 1
      nyA         = 3
      nA          = nxA + nyA + 1      # No of free parameters/Conditions to satisfy
      
      FAC         = 1.0

      # Incident points     
      mA          = 3       
      xA          = zeros(Float64,mA)
      yA          = zeros(Float64,mA)
      if mA > 0
        xA        = FAC*[0.0; -0.15;   2.0]
        yA        = FAC*[0.0;  0.32;   4.0]
      end  
      
      # X-Derivative Points
      mxA = 0
      xdxA = zeros(Float64,mxA)
      ydxA = zeros(Float64,mxA)
      if (mxA>0)
        xdxA[1] = xA[2]
        ydxA[1] = yA[2] 
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
      #-------------------- 
      nxB               = 1
      nyB               = 1
      nB                = nxB + nyB + 1      # No of free parameters/Conditions to satisfy
      
      # Incident points
      mB                = 2
      xB                = zeros(Float64,mB)
      yB                = zeros(Float64,mB)
      if mB > 0
        xB              = FAC*[0.0;  3.0]
        yB              = FAC*[0.0;  4.5]
      end 
      
      # X-Derivative Points
      mxB = 0
      xdxB              = zeros(Float64,mxB)
      ydxB              = zeros(Float64,mxB)
      if (mxB>0)
        xdxB[1]            = xB[2]
        ydxB[1]            = yB[2] 
      end  
      
      # Y-Derivative Points
      myB = 0
      xdyB              = zeros(Float64,myB)
      ydyB              = zeros(Float64,myB)
      
      if myB>0
        xdyB[1]         = xB[2]
        ydyB[1]         = yB[2]
      end  

    elseif nullcline_set == 23
#----------------------------------------       
      println("Testing Turing Patterns")

      nxA         = 1
      nyA         = 3
      nA          = nxA + nyA + 1      # No of free parameters/Conditions to satisfy
      
      FAC         = 1.0

      # Incident points     
      mA          = 3       
      xA          = zeros(Float64,mA)
      yA          = zeros(Float64,mA)
      if mA > 0
        xA        = FAC*[0.0; -0.01;   1.0]
        yA        = FAC*[0.0;  0.02;   2.0]
      end  
      
      # X-Derivative Points
      mxA = 0
      xdxA = zeros(Float64,mxA)
      ydxA = zeros(Float64,mxA)
      if (mxA>0)
        xdxA[1] = xA[2]
        ydxA[1] = yA[2] 
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
      #-------------------- 
      nxB               = 1
      nyB               = 1
      nB                = nxB + nyB + 1      # No of free parameters/Conditions to satisfy
      
      # Incident points
      mB                = 2
      xB                = zeros(Float64,mB)
      yB                = zeros(Float64,mB)
      if mB > 0
        xB              = FAC*[0.0;  2.5]
        yB              = FAC*[0.0;  0.05]
      end 
      
      # X-Derivative Points
      mxB = 0
      xdxB              = zeros(Float64,mxB)
      ydxB              = zeros(Float64,mxB)
      if (mxB>0)
        xdxB[1]            = xB[2]
        ydxB[1]            = yB[2] 
      end  
      
      # Y-Derivative Points
      myB = 0
      xdyB              = zeros(Float64,myB)
      ydyB              = zeros(Float64,myB)
      
      if myB>0
        xdyB[1]         = xB[2]
        ydyB[1]         = yB[2]
      end  

    elseif nullcline_set == 24
#----------------------------------------       
      println("Symmetric fixed points on the upper and lower branches")

      nxA         = 1
      nyA         = 3
      nA          = nxA + nyA + 1      # No of free parameters/Conditions to satisfy
      
      FAC         = 1.0

      # Incident points     
      mA          = 3       
      xA          = zeros(Float64,mA)
      yA          = zeros(Float64,mA)
      if mA > 0
        xA        = FAC*[ 0.0;  -3.0;  3.0]
        yA        = FAC*[-1.0;  -4.0;  2.0]
      end
      
      # X-Derivative Points
      mxA = 0
      xdxA = zeros(Float64,mxA)
      ydxA = zeros(Float64,mxA)
      if (mxA>0)
        xdxA[1] = xA[2]
        ydxA[1] = yA[2] 
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
      #-------------------- 
      nxB               = 1
      nyB               = 1
      nB                = nxB + nyB + 1      # No of free parameters/Conditions to satisfy
      
      # Incident points
      mB                = 2
      xB                = zeros(Float64,mB)
      yB                = zeros(Float64,mB)
      if mB > 0
        xB              = FAC*[0.0;  1.0]
        yB              = FAC*[0.0; -8.0]
       
      end  
      
      # X-Derivative Points
      mxB = 0
      xdxB              = zeros(Float64,mxB)
      ydxB              = zeros(Float64,mxB)
      if (mxB>0)
        xdxB[1]            = xB[2]
        ydxB[1]            = yB[2] 
      end  
      
      # Y-Derivative Points
      myB = 0
      xdyB              = zeros(Float64,myB)
      ydyB              = zeros(Float64,myB)
      
      if myB>0
        xdyB[1]         = xB[2]
        ydyB[1]         = yB[2]
      end

    elseif nullcline_set == 30
#----------------------------------------       
      println("Slugs. Small non-linear instability threshold")

      nxA         = 1
      nyA         = 3
      nA          = nxA + nyA + 1      # No of free parameters/Conditions to satisfy
      
      FAC         = 1.0

      # Incident points     
      mA          = 3       
      xA          = zeros(Float64,mA)
      yA          = zeros(Float64,mA)
      if mA > 0
        xA        = FAC*[0.0; -0.15;   2.0]
        yA        = FAC*[0.0;  0.32;   4.0]
      end  
      
      # X-Derivative Points
      mxA = 0
      xdxA = zeros(Float64,mxA)
      ydxA = zeros(Float64,mxA)
      if (mxA>0)
        xdxA[1] = xA[2]
        ydxA[1] = yA[2] 
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
      #-------------------- 
      nxB               = 1
      nyB               = 1
      nB                = nxB + nyB + 1      # No of free parameters/Conditions to satisfy
      
      # Incident points
      mB                = 2
      xB                = zeros(Float64,mB)
      yB                = zeros(Float64,mB)
      if mB > 0
        xB              = FAC*[0.0;  3.0]
        yB              = FAC*[0.0;  4.285]
      end 
      
      # X-Derivative Points
      mxB = 0
      xdxB              = zeros(Float64,mxB)
      ydxB              = zeros(Float64,mxB)
      if (mxB>0)
        xdxB[1]            = xB[2]
        ydxB[1]            = yB[2] 
      end  
      
      # Y-Derivative Points
      myB = 0
      xdyB              = zeros(Float64,myB)
      ydyB              = zeros(Float64,myB)
      
      if myB>0
        xdyB[1]         = xB[2]
        ydyB[1]         = yB[2]
      end  


    elseif nullcline_set == 51
#----------------------------------------       
      println("Nullcline for dynamic switching")

      nxA         = 1
      nyA         = 3
      nA          = nxA + nyA + 1      # No of free parameters/Conditions to satisfy
      
      FAC         = 1.0

      # Incident points     
      mA          = 3 
      xA          = zeros(Float64,mA)
      yA          = zeros(Float64,mA)
      if mA > 0
        xA        = FAC*[1.0;   6.0;  -4.0]
        yA        = FAC*[0.0;   0.7;  -0.7]
      end  
      
      # X-Derivative Points
      mxA = 0
      xdxA = zeros(Float64,mxA)
      ydxA = zeros(Float64,mxA)
      if (mxA>0)
        xdxA[1] = xA[3]
        ydxA[1] = yA[3] 
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
      #-------------------- 
      nxB               = 1
      nyB               = 2
      nB                = nxB + nyB + 1      # No of free parameters/Conditions to satisfy
      
      # Incident points
      mB                = 2
      xB                = zeros(Float64,mB)
      yB                = zeros(Float64,mB)
      if mB > 0
        xB              = FAC*[0.0;  1.0]
        yB              = FAC*[0.0;  3.0]
      end  
      
      # X-Derivative Points
      mxB = 0
      xdxB              = zeros(Float64,mxB)
      ydxB              = zeros(Float64,mxB)
      if (mxB>0)
        xdxB[1]            = xB[2]
        ydxB[1]            = yB[2] 
      end  
      
      # Y-Derivative Points
      myB = 1
      xdyB              = zeros(Float64,myB)
      ydyB              = zeros(Float64,myB)
      if myB>0
        xdyB[1]         = xB[2]
        ydyB[1]         = yB[2]
      end  
    elseif nullcline_set == 52
#----------------------------------------       
      println("Nullcline for dynamic switching 2")

      nxA         = 1
      nyA         = 3
      nA          = nxA + nyA + 1      # No of free parameters/Conditions to satisfy
      
      FAC         = 1.0

      # Incident points     
      mA          = 3 
      xA          = zeros(Float64,mA)
      yA          = zeros(Float64,mA)
      if mA > 0
        xA        = FAC*[1.0;     6.0]
        yA        = FAC*[0.165;    0.835]
      end  
      
      # X-Derivative Points
      mxA = 0
      xdxA = zeros(Float64,mxA)
      ydxA = zeros(Float64,mxA)
      if (mxA>0)
        xdxA[1] = xA[3]
        ydxA[1] = yA[3] 
      end  
      
      # Y-Derivative Points
      myA = 2
      xdyA = zeros(Float64,myA)
      ydyA = zeros(Float64,myA)
      if myA>0
        xdyA[1] = xA[1]
        ydyA[1] = yA[1]

        xdyA[2] = xA[2]
        ydyA[2] = yA[2]
      end  
      
      # Points for G(x,y)
      #-------------------- 
      nxB               = 1
      nyB               = 2
      nB                = nxB + nyB + 1      # No of free parameters/Conditions to satisfy
      
      # Incident points
      mB                = 2
      xB                = zeros(Float64,mB)
      yB                = zeros(Float64,mB)
      if mB > 0
        xB              = FAC*[0.0;  1.0]
        yB              = FAC*[0.0;  3.0]
      end  
      
      # X-Derivative Points
      mxB = 0
      xdxB              = zeros(Float64,mxB)
      ydxB              = zeros(Float64,mxB)
      if (mxB>0)
        xdxB[1]            = xB[2]
        ydxB[1]            = yB[2] 
      end  
      
      # Y-Derivative Points
      myB = 1
      xdyB              = zeros(Float64,myB)
      ydyB              = zeros(Float64,myB)
      if myB>0
        xdyB[1]         = xB[2]
        ydyB[1]         = yB[2]
      end  
    elseif nullcline_set == 53
#----------------------------------------       
      println("Nullcline for dynamic switching 3")

      nxA         = 1
      nyA         = 3
      nA          = nxA + nyA + 1      # No of free parameters/Conditions to satisfy
      
      FAC         = 1.0

      # Incident points     
      mA          = 2 
      xA          = zeros(Float64,mA)
      yA          = zeros(Float64,mA)
      if mA > 0
        xA        = FAC*[1.60;  5.00]
        yA        = FAC*[0.21;  0.75]
      end  
      
      # X-Derivative Points
      mxA = 0 
      xdxA = zeros(Float64,mxA)
      ydxA = zeros(Float64,mxA)
      if (mxA>0)
        xdxA[1] = xA[3]
        ydxA[1] = yA[3] 
      end  
      
      # Y-Derivative Points
      myA = 2
      xdyA = zeros(Float64,myA)
      ydyA = zeros(Float64,myA)
      if myA>0
        xdyA[1] = xA[1]
        ydyA[1] = yA[1]

        xdyA[2] = xA[2]
        ydyA[2] = yA[2]
      end  
      
      # Points for G(x,y)
      #-------------------- 
      nxB               = 1
      nyB               = 2
      nB                = nxB + nyB + 1      # No of free parameters/Conditions to satisfy
      
      # Incident points
      mB                = 2
      xB                = zeros(Float64,mB)
      yB                = zeros(Float64,mB)
      if mB > 0
        xB              = FAC*[0.0;  1.0]
        yB              = FAC*[0.0;  3.0]
      end  
      
      # X-Derivative Points
      mxB = 0
      xdxB              = zeros(Float64,mxB)
      ydxB              = zeros(Float64,mxB)
      if (mxB>0)
        xdxB[1]            = xB[2]
        ydxB[1]            = yB[2] 
      end  
      
      # Y-Derivative Points
      myB = 1
      xdyB              = zeros(Float64,myB)
      ydyB              = zeros(Float64,myB)
      if myB>0
        xdyB[1]         = xB[2]
        ydyB[1]         = yB[2]
      end  
    elseif nullcline_set == 54
#----------------------------------------       
      println("Nullcline for dynamic switching 4")

      nxA         = 1
      nyA         = 3
      nA          = nxA + nyA + 1      # No of free parameters/Conditions to satisfy
      
      FAC         = 1.0

      # Incident points     
      mA          = 3 
      xA          = zeros(Float64,mA)
      yA          = zeros(Float64,mA)
      if mA > 0 
        xA        = [ 0.20;  0.40;  0.30]
        yA        = [-0.90;  0.90;  0.00]
      end  
      
      # X-Derivative Points
      mxA = 0 
      xdxA = zeros(Float64,mxA)
      ydxA = zeros(Float64,mxA)
      if (mxA>0)
        xdxA[1] = xA[2]
        ydxA[1] = yA[2] 
      end  
      
      # Y-Derivative Points
      myA = 1
      xdyA = zeros(Float64,myA)
      ydyA = zeros(Float64,myA)
      if myA>0
        xdyA[1] = xA[1]
        ydyA[1] = yA[1]

        # xdyA[2] = xA[2]
        # ydyA[2] = yA[2]
      end  
      
      # Points for G(x,y)
      #-------------------- 
      nxB               = 1
      nyB               = 2
      nB                = nxB + nyB + 1      # No of free parameters/Conditions to satisfy
      
      # Incident points
      mB                = 2
      xB                = zeros(Float64,mB)
      yB                = zeros(Float64,mB)
      if mB > 0
        xB              = FAC*[0.0;  1.0]
        yB              = FAC*[0.0;  3.0]
      end  
      
      # X-Derivative Points
      mxB = 0
      xdxB              = zeros(Float64,mxB)
      ydxB              = zeros(Float64,mxB)
      if (mxB>0)
        xdxB[1]            = xB[2]
        ydxB[1]            = yB[2] 
      end  
      
      # Y-Derivative Points
      myB = 1
      xdyB              = zeros(Float64,myB)
      ydyB              = zeros(Float64,myB)
      if myB>0
        xdyB[1]         = xB[2]
        ydyB[1]         = yB[2]
      end  

    elseif nullcline_set == 55
#----------------------------------------       
      println("Nullcline for dynamic switching 5")

      nxA         = 1
      nyA         = 3
      nA          = nxA + nyA + 1      # No of free parameters/Conditions to satisfy
      
      FAC         = 1.0

      # Incident points     
      mA          = 3 
      xA          = zeros(Float64,mA)
      yA          = zeros(Float64,mA)
      if mA > 0 
        xA        = [ 0.40;  0.20;  0.30]
        yA        = [-0.90;  0.90;  0.00]
      end  
      
      # X-Derivative Points
      mxA = 0 
      xdxA = zeros(Float64,mxA)
      ydxA = zeros(Float64,mxA)
      if (mxA>0)
        xdxA[1] = xA[2]
        ydxA[1] = yA[2] 
      end  
      
      # Y-Derivative Points
      myA = 1
      xdyA = zeros(Float64,myA)
      ydyA = zeros(Float64,myA)
      if myA>0
        xdyA[1] = xA[1]
        ydyA[1] = yA[1]

        # xdyA[2] = xA[2]
        # ydyA[2] = yA[2]
      end  
      
      # Points for G(x,y)
      #-------------------- 
      nxB               = 1
      nyB               = 2
      nB                = nxB + nyB + 1      # No of free parameters/Conditions to satisfy
      
      # Incident points
      mB                = 2
      xB                = zeros(Float64,mB)
      yB                = zeros(Float64,mB)
      if mB > 0
        xB              = FAC*[0.0;  1.0]
        yB              = FAC*[0.0;  3.0]
      end  
      
      # X-Derivative Points
      mxB = 0
      xdxB              = zeros(Float64,mxB)
      ydxB              = zeros(Float64,mxB)
      if (mxB>0)
        xdxB[1]            = xB[2]
        ydxB[1]            = yB[2] 
      end  
      
      # Y-Derivative Points
      myB = 1
      xdyB              = zeros(Float64,myB)
      ydyB              = zeros(Float64,myB)
      if myB>0
        xdyB[1]         = xB[2]
        ydyB[1]         = yB[2]
      end  

    elseif nullcline_set == 56
#----------------------------------------       
      println("Nullcline for dynamic switching 6: Origin")

      nxA         = 1
      nyA         = 3
      nA          = nxA + nyA + 1      # No of free parameters/Conditions to satisfy
      
      FAC         = 1.0

      # Incident points     
      mA          = 3 
      xA          = zeros(Float64,mA)
      yA          = zeros(Float64,mA)
      if mA > 0 
        xA        = [ 0.10; -0.10;  0.00]
        yA        = [-0.90;  0.90;  0.00]
      end  
      
      # X-Derivative Points
      mxA = 0 
      xdxA = zeros(Float64,mxA)
      ydxA = zeros(Float64,mxA)
      if (mxA>0)
        xdxA[1] = xA[2]
        ydxA[1] = yA[2] 
      end  
      
      # Y-Derivative Points
      myA = 1
      xdyA = zeros(Float64,myA)
      ydyA = zeros(Float64,myA)
      if myA>0
        xdyA[1] = xA[1]
        ydyA[1] = yA[1]

        # xdyA[2] = xA[2]
        # ydyA[2] = yA[2]
      end  
      
      # Points for G(x,y)
      #-------------------- 
      nxB               = 1
      nyB               = 2
      nB                = nxB + nyB + 1      # No of free parameters/Conditions to satisfy
      
      # Incident points
      mB                = 2
      xB                = zeros(Float64,mB)
      yB                = zeros(Float64,mB)
      if mB > 0
        xB              = FAC*[0.0;  1.0]
        yB              = FAC*[0.0;  3.0]
      end  
      
      # X-Derivative Points
      mxB = 0
      xdxB              = zeros(Float64,mxB)
      ydxB              = zeros(Float64,mxB)
      if (mxB>0)
        xdxB[1]            = xB[2]
        ydxB[1]            = yB[2] 
      end  
      
      # Y-Derivative Points
      myB = 1
      xdyB              = zeros(Float64,myB)
      ydyB              = zeros(Float64,myB)
      if myB>0
        xdyB[1]         = xB[2]
        ydyB[1]         = yB[2]
      end  

    elseif nullcline_set == 101
#----------------------------------------       
      println("Slugs with large threshold ")

      nxA         = 1
      nyA         = 3
      nA          = nxA + nyA + 1      # No of free parameters/Conditions to satisfy
      
      FAC         = 1.0

      # Incident points     
      mA          = 3       
      xA          = zeros(Float64,mA)
      yA          = zeros(Float64,mA)
      if mA > 0
        xA        = FAC*[0.0; -1.0;    3.0]
        yA        = FAC*[0.0;  0.8;    5.0]
      end  
      
      # X-Derivative Points
      mxA = 0
      xdxA = zeros(Float64,mxA)
      ydxA = zeros(Float64,mxA)
      if (mxA>0)
        xdxA[1] = xA[2]
        ydxA[1] = yA[2] 
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
      #-------------------- 
      nxB               = 1
      nyB               = 2
      nB                = nxB + nyB + 1      # No of free parameters/Conditions to satisfy
      
      # Incident points
      mB                = 2
      xB                = zeros(Float64,mB)
      yB                = zeros(Float64,mB)
      if mB > 0
        xB              = FAC*[0.0;  2.5]
        yB              = FAC*[0.0;  5.0]
      end  
      
      # X-Derivative Points
      mxB = 0
      xdxB              = zeros(Float64,mxB)
      ydxB              = zeros(Float64,mxB)
      if (mxB>0)
        xdxB[1]            = xB[2]
        ydxB[1]            = yB[2] 
      end  
      
      # Y-Derivative Points
      myB = 1
      xdyB              = zeros(Float64,myB)
      ydyB              = zeros(Float64,myB)
      
      if myB>0
        xdyB[1]         = xB[2]
        ydyB[1]         = yB[2]
      end  

#     Paper parameters
#--------------------------------------------------     

    elseif nullcline_set == 201
#----------------------------------------       
      println("Paper Figure 1")
      println("Paper Wedges")

      nxA         = 1
      nyA         = 3
      nA          = nxA + nyA + 1      # No of free parameters/Conditions to satisfy
      
      FAC         = 1.0

      # Incident points     
      mA          = 3       
      xA          = zeros(Float64,mA)
      yA          = zeros(Float64,mA)
      if mA > 0
        xA        = FAC*[0.0; -0.15;   2.0]
        yA        = FAC*[0.0;  0.32;   4.0]
      end  
      
      # X-Derivative Points
      mxA = 0
      xdxA = zeros(Float64,mxA)
      ydxA = zeros(Float64,mxA)
      if (mxA>0)
        xdxA[1] = xA[2]
        ydxA[1] = yA[2] 
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
      #-------------------- 
      nxB               = 1
      nyB               = 1
      nB                = nxB + nyB + 1      # No of free parameters/Conditions to satisfy
      
      # Incident points
      mB                = 2
      xB                = zeros(Float64,mB)
      yB                = zeros(Float64,mB)
      if mB > 0
        xB              = FAC*[0.0;  5.0]
        yB              = FAC*[0.0;  4.0]
      end 
      
      # X-Derivative Points
      mxB = 0
      xdxB              = zeros(Float64,mxB)
      ydxB              = zeros(Float64,mxB)
      if (mxB>0)
        xdxB[1]            = xB[2]
        ydxB[1]            = yB[2] 
      end  
      
      # Y-Derivative Points
      myB = 0
      xdyB              = zeros(Float64,myB)
      ydyB              = zeros(Float64,myB)
      
      if myB>0
        xdyB[1]         = xB[2]
        ydyB[1]         = yB[2]
      end  

    elseif nullcline_set == 202
#----------------------------------------       
      println("Symmetric Limit Cycle Oscillations (LCOs)")

      nxA         = 1
      nyA         = 3
      nA          = nxA + nyA + 1      # No of free parameters/Conditions to satisfy
      
      FAC         = 1.0

      # Incident points     
      mA          = 3       
      xA          = zeros(Float64,mA)
      yA          = zeros(Float64,mA)
      if mA > 0
        xA        = FAC*[0.0; -0.15;   2.0]
        yA        = FAC*[0.0;  0.32;   4.0]
      end  
      
      # X-Derivative Points
      mxA = 0
      xdxA = zeros(Float64,mxA)
      ydxA = zeros(Float64,mxA)
      if (mxA>0)
        xdxA[1] = xA[2]
        ydxA[1] = yA[2] 
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
      #-------------------- 
      nxB               = 1
      nyB               = 1
      nB                = nxB + nyB + 1      # No of free parameters/Conditions to satisfy
      
      # Incident points
      mB                = 2
      xB                = zeros(Float64,mB)
      yB                = zeros(Float64,mB)
      if mB > 0
        xB              = FAC*[1.6;    4.5]
        yB              = FAC*[1.71;  2.70]
      end 
      
      # X-Derivative Points
      mxB = 0
      xdxB              = zeros(Float64,mxB)
      ydxB              = zeros(Float64,mxB)
      if (mxB>0)
        xdxB[1]            = xB[2]
        ydxB[1]            = yB[2] 
      end  
      
      # Y-Derivative Points
      myB = 0
      xdyB              = zeros(Float64,myB)
      ydyB              = zeros(Float64,myB)
      
      if myB>0
        xdyB[1]         = xB[2]
        ydyB[1]         = yB[2]
      end  

    elseif nullcline_set == 203
#----------------------------------------       
      println("Upper-branch dominated Limit Cycle Oscillations (LCOs)")

      nxA         = 1
      nyA         = 3
      nA          = nxA + nyA + 1      # No of free parameters/Conditions to satisfy
      
      FAC         = 1.0

      # Incident points     
      mA          = 3       
      xA          = zeros(Float64,mA)
      yA          = zeros(Float64,mA)
      if mA > 0
        xA        = FAC*[0.0; -0.15;   2.0]
        yA        = FAC*[0.0;  0.32;   4.0]
      end  
      
      # X-Derivative Points
      mxA = 0
      xdxA = zeros(Float64,mxA)
      ydxA = zeros(Float64,mxA)
      if (mxA>0)
        xdxA[1] = xA[2]
        ydxA[1] = yA[2] 
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
      #-------------------- 
      nxB               = 1
      nyB               = 1
      nB                = nxB + nyB + 1      # No of free parameters/Conditions to satisfy
      
      # Incident points
      mB                = 2
      xB                = zeros(Float64,mB)
      yB                = zeros(Float64,mB)
      if mB > 0
        xB              = FAC*[1.6;    4.5]
        yB              = FAC*[2.01;  3.00]
      end 
      
      # X-Derivative Points
      mxB = 0
      xdxB              = zeros(Float64,mxB)
      ydxB              = zeros(Float64,mxB)
      if (mxB>0)
        xdxB[1]            = xB[2]
        ydxB[1]            = yB[2] 
      end  
      
      # Y-Derivative Points
      myB = 0
      xdyB              = zeros(Float64,myB)
      ydyB              = zeros(Float64,myB)
      
      if myB>0
        xdyB[1]         = xB[2]
        ydyB[1]         = yB[2]
      end  

    elseif nullcline_set == 204
#----------------------------------------       
      println("Lower-branch dominated Limit Cycle Oscillations (LCOs)")

      nxA         = 1
      nyA         = 3
      nA          = nxA + nyA + 1      # No of free parameters/Conditions to satisfy
      
      FAC         = 1.0

      # Incident points     
      mA          = 3       
      xA          = zeros(Float64,mA)
      yA          = zeros(Float64,mA)
      if mA > 0
        xA        = FAC*[0.0; -0.15;   2.0]
        yA        = FAC*[0.0;  0.32;   4.0]
      end  
      
      # X-Derivative Points
      mxA = 0
      xdxA = zeros(Float64,mxA)
      ydxA = zeros(Float64,mxA)
      if (mxA>0)
        xdxA[1] = xA[2]
        ydxA[1] = yA[2] 
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
      #-------------------- 
      nxB               = 1
      nyB               = 1
      nB                = nxB + nyB + 1      # No of free parameters/Conditions to satisfy
      
      # Incident points
      mB                = 2
      xB                = zeros(Float64,mB)
      yB                = zeros(Float64,mB)
      if mB > 0
        xB              = FAC*[1.6;    4.5]
        yB              = FAC*[0.96;  1.95]
      end 
      
      # X-Derivative Points
      mxB = 0
      xdxB              = zeros(Float64,mxB)
      ydxB              = zeros(Float64,mxB)
      if (mxB>0)
        xdxB[1]            = xB[2]
        ydxB[1]            = yB[2] 
      end  
      
      # Y-Derivative Points
      myB = 0
      xdyB              = zeros(Float64,myB)
      ydyB              = zeros(Float64,myB)
      
      if myB>0
        xdyB[1]         = xB[2]
        ydyB[1]         = yB[2]
      end  

    elseif nullcline_set == 205
#----------------------------------------       
      println("Slugs.")

      nxA         = 1
      nyA         = 3
      nA          = nxA + nyA + 1      # No of free parameters/Conditions to satisfy
      
      FAC         = 1.0

      # Incident points     
      mA          = 3       
      xA          = zeros(Float64,mA)
      yA          = zeros(Float64,mA)
      if mA > 0
        xA        = FAC*[0.0; -0.15;   2.0]
        yA        = FAC*[0.0;  0.32;   4.0]
      end  
      
      # X-Derivative Points
      mxA = 0
      xdxA = zeros(Float64,mxA)
      ydxA = zeros(Float64,mxA)
      if (mxA>0)
        xdxA[1] = xA[2]
        ydxA[1] = yA[2] 
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
      #-------------------- 
      nxB               = 1
      nyB               = 1
      nB                = nxB + nyB + 1      # No of free parameters/Conditions to satisfy
      
      # Incident points
      mB                = 2
      xB                = zeros(Float64,mB)
      yB                = zeros(Float64,mB)
      if mB > 0
        xB              = FAC*[0.0;  2.5]
        yB              = FAC*[0.0;  4.0]
      end 
      
      # X-Derivative Points
      mxB = 0
      xdxB              = zeros(Float64,mxB)
      ydxB              = zeros(Float64,mxB)
      if (mxB>0)
        xdxB[1]            = xB[2]
        ydxB[1]            = yB[2] 
      end  
      
      # Y-Derivative Points
      myB = 0
      xdyB              = zeros(Float64,myB)
      ydyB              = zeros(Float64,myB)
      
      if myB>0
        xdyB[1]         = xB[2]
        ydyB[1]         = yB[2]
      end  

    elseif nullcline_set == 206
#----------------------------------------       
      println("Stripes")

      nxA         = 1
      nyA         = 3
      nA          = nxA + nyA + 1      # No of free parameters/Conditions to satisfy
      
      FAC         = 1.0

      # Incident points     
      mA          = 3       
      xA          = zeros(Float64,mA)
      yA          = zeros(Float64,mA)
      if mA > 0
        xA        = FAC*[0.0; -0.15;   2.0]
        yA        = FAC*[0.0;  0.32;   4.0]
      end  
      
      # X-Derivative Points
      mxA = 0
      xdxA = zeros(Float64,mxA)
      ydxA = zeros(Float64,mxA)
      if (mxA>0)
        xdxA[1] = xA[2]
        ydxA[1] = yA[2] 
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
      #-------------------- 
      nxB               = 1
      nyB               = 1
      nB                = nxB + nyB + 1      # No of free parameters/Conditions to satisfy
      
      # Incident points
      mB                = 2
      xB                = zeros(Float64,mB)
      yB                = zeros(Float64,mB)
      if mB > 0
        xB              = FAC*[1.6;    2.5]
        yB              = FAC*[1.71;   0.0]
      end 
      
      # X-Derivative Points
      mxB = 0
      xdxB              = zeros(Float64,mxB)
      ydxB              = zeros(Float64,mxB)
      if (mxB>0)
        xdxB[1]            = xB[2]
        ydxB[1]            = yB[2] 
      end  
      
      # Y-Derivative Points
      myB = 0
      xdyB              = zeros(Float64,myB)
      ydyB              = zeros(Float64,myB)
      
      if myB>0
        xdyB[1]         = xB[2]
        ydyB[1]         = yB[2]
      end  

    elseif nullcline_set == 207
#----------------------------------------       
      println("Branching")

      nxA         = 1
      nyA         = 3
      nA          = nxA + nyA + 1      # No of free parameters/Conditions to satisfy
      
      FAC         = 1.0

      # Incident points     
      mA          = 3       
      xA          = zeros(Float64,mA)
      yA          = zeros(Float64,mA)
      if mA > 0
        xA        = FAC*[0.0; -0.15;   2.0]
        yA        = FAC*[0.0;  0.32;   4.0]
      end  
      
      # X-Derivative Points
      mxA = 0
      xdxA = zeros(Float64,mxA)
      ydxA = zeros(Float64,mxA)
      if (mxA>0)
        xdxA[1] = xA[2]
        ydxA[1] = yA[2] 
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
      #-------------------- 
      nxB               = 1
      nyB               = 1
      nB                = nxB + nyB + 1      # No of free parameters/Conditions to satisfy
      
      # Incident points
      mB                = 2
      xB                = zeros(Float64,mB)
      yB                = zeros(Float64,mB)
      if mB > 0
#        xB              = FAC*[0.0;  4.0]
#        yB              = FAC*[0.0;  4.0]

        xB              = FAC*[0.0;  2.5]
        yB              = FAC*[0.0;  4.25]
      end 
      
      # X-Derivative Points
      mxB = 0
      xdxB              = zeros(Float64,mxB)
      ydxB              = zeros(Float64,mxB)
      if (mxB>0)
        xdxB[1]            = xB[2]
        ydxB[1]            = yB[2] 
      end  
      
      # Y-Derivative Points
      myB = 0
      xdyB              = zeros(Float64,myB)
      ydyB              = zeros(Float64,myB)
      
      if myB>0
        xdyB[1]         = xB[2]
        ydyB[1]         = yB[2]
      end  

    elseif nullcline_set == 208
#----------------------------------------       
      println("Branching/Crossings")

      nxA         = 1
      nyA         = 3
      nA          = nxA + nyA + 1      # No of free parameters/Conditions to satisfy
      
      FAC         = 1.0

      # Incident points     
      mA          = 3       
      xA          = zeros(Float64,mA)
      yA          = zeros(Float64,mA)
      if mA > 0
        xA        = FAC*[0.0; -0.15;   2.0]
        yA        = FAC*[0.0;  0.32;   4.0]
      end  
      
      # X-Derivative Points
      mxA = 0
      xdxA = zeros(Float64,mxA)
      ydxA = zeros(Float64,mxA)
      if (mxA>0)
        xdxA[1] = xA[2]
        ydxA[1] = yA[2] 
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
      #-------------------- 
      nxB               = 1
      nyB               = 1
      nB                = nxB + nyB + 1      # No of free parameters/Conditions to satisfy
      
      # Incident points
      mB                = 2
      xB                = zeros(Float64,mB)
      yB                = zeros(Float64,mB)
      if mB > 0
        xB              = FAC*[0.0;  3.00]
        yB              = FAC*[0.0;  4.285]

#        xB              = FAC*[0.0;  4.25]
#        yB              = FAC*[0.0;  5.0]
      end 
      
      # X-Derivative Points
      mxB = 0
      xdxB              = zeros(Float64,mxB)
      ydxB              = zeros(Float64,mxB)
      if (mxB>0)
        xdxB[1]            = xB[2]
        ydxB[1]            = yB[2] 
      end  
      
      # Y-Derivative Points
      myB = 0
      xdyB              = zeros(Float64,myB)
      ydyB              = zeros(Float64,myB)
      
      if myB>0
        xdyB[1]         = xB[2]
        ydyB[1]         = yB[2]
      end  

    elseif nullcline_set == 209
#----------------------------------------       
      println("One Sided Triangles")

      nxA         = 1
      nyA         = 3
      nA          = nxA + nyA + 1      # No of free parameters/Conditions to satisfy
      
      FAC         = 1.0

      # Incident points     
      mA          = 3       
      xA          = zeros(Float64,mA)
      yA          = zeros(Float64,mA)
      if mA > 0
#        xA        = FAC*[0.0; -0.15;   2.0]
#        yA        = FAC*[0.0;  0.32;   4.0]

        xA        = FAC*[0.0; -0.5;    2.0]
        yA        = FAC*[0.0;  0.6;    5.0]
      end  
      
      # X-Derivative Points
      mxA = 0
      xdxA = zeros(Float64,mxA)
      ydxA = zeros(Float64,mxA)
      if (mxA>0)
        xdxA[1] = xA[2]
        ydxA[1] = yA[2] 
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
      #-------------------- 
      nxB               = 1
      nyB               = 1
      nB                = nxB + nyB + 1      # No of free parameters/Conditions to satisfy
      
      # Incident points
      mB                = 2
      xB                = zeros(Float64,mB)
      yB                = zeros(Float64,mB)
      if mB > 0
        xB              = FAC*[0.0;  3.0]
        yB              = FAC*[0.0;  4.75]

#        xB              = FAC*[0.0;  4.1]
#        yB              = FAC*[0.0;  4.6]
      end 
      
      # X-Derivative Points
      mxB = 0
      xdxB              = zeros(Float64,mxB)
      ydxB              = zeros(Float64,mxB)
      if (mxB>0)
        xdxB[1]            = xB[2]
        ydxB[1]            = yB[2] 
      end  
      
      # Y-Derivative Points
      myB = 0
      xdyB              = zeros(Float64,myB)
      ydyB              = zeros(Float64,myB)
      
      if myB>0
        xdyB[1]         = xB[2]
        ydyB[1]         = yB[2]
      end  
    elseif nullcline_set == 210
#----------------------------------------       
      println("Sierpinsky Triangles")

      nxA         = 1
      nyA         = 3
      nA          = nxA + nyA + 1      # No of free parameters/Conditions to satisfy
      
      FAC         = 1.0

      # Incident points     
      mA          = 3       
      xA          = zeros(Float64,mA)
      yA          = zeros(Float64,mA)
      if mA > 0
#        xA        = FAC*[0.0; -0.15;   2.0]
#        yA        = FAC*[0.0;  0.32;   4.0]

        xA        = FAC*[0.0; -0.5;    3.5]
        yA        = FAC*[0.0;  0.6;    4.6]
      end  
      
      # X-Derivative Points
      mxA = 0
      xdxA = zeros(Float64,mxA)
      ydxA = zeros(Float64,mxA)
      if (mxA>0)
        xdxA[1] = xA[2]
        ydxA[1] = yA[2] 
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
      #-------------------- 
      nxB               = 1
      nyB               = 1
      nB                = nxB + nyB + 1      # No of free parameters/Conditions to satisfy
      
      # Incident points
      mB                = 2
      xB                = zeros(Float64,mB)
      yB                = zeros(Float64,mB)
      if mB > 0
        xB              = FAC*[0.0;  3.5]
        yB              = FAC*[0.0;  4.2]
      end 
      
      # X-Derivative Points
      mxB = 0
      xdxB              = zeros(Float64,mxB)
      ydxB              = zeros(Float64,mxB)
      if (mxB>0)
        xdxB[1]            = xB[2]
        ydxB[1]            = yB[2] 
      end  
      
      # Y-Derivative Points
      myB = 0
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
@doc """

   function GetNullClineParams(nullcline_set)

            Null-cline Sets: 1         : Pulses
                           : 2         : Slugs. No change in instability threshold
                           : 3         : Pulses. Smaller instability threshold
                           : 4         : Slugs. Smaller instability threshold
                           : 5         : Extreme Slugs. Smaller instability threshold
                           : 6         : Extreme Pulse collapse
                           : 7         : Extreme slugs
                           : 8         : Weak slugs
                           : 9         : Slugs. Linear Inhibitor.
                           : 10        : Limit-cycline Oscillation. Activation dominated
                           : 11        : Limit-cycline Oscillation. De-activation dominated
                           : 12        : Limit-cycline Oscillation. De-activation dominated (cubic in x)
                           : 13        : Two fixed points - upper and lower branch.
                           : 14        : Symmetric fixed points - upper and lower branch
                           : 15        : Symmetric LCO
                           : 16        : Two Asymmetric fixed points
                           : 17        : Unstable G-null-cline
                           : 18        : Marginal Instability
                           : 19        : Stable Oscillatory Linear mode
                           : 20        : Slugs. Large non-linear instability threshold
                           : 21        : Slugs. Small non-linear instability threshold
                           : 51        : Dynamic Switching (λ)
                           : 52        : Dynamic Switching (λ)
                           : 53        : Dynamic Switching (λ)
                           : 54        : Dynamic Switching (λ)
                           : 55        : Dynamic Switching (λ)
                           :
                           :101        : Slug with large threshold

               Output the parameters as a NullClineParam Structure.                  

"""
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
    
    return params 
end
#---------------------------------------------------------------------- 


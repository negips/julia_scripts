#!/bin/julia
# Build the Evolution functions based on the parameters
#---------------------------------------------------------------------- 
"""
   function FX(x,c0,cx)

     F(x) = c0 + Σcx_i*x^i

"""
function FX(x,c0,cx)
  nx = length(cx)

  s = c0
  for i in 1:nx
    s = s .+ cx[i]*x.^i
  end

  return s
end  
#---------------------------------------------------------------------- 
"""
   function GradFX(x,cx)

     ∇(F(x)) = Σ(cx_i)*i*x^(i-1)

"""
function GradFX(x,cx)
  nx = length(cx)

  s = 0.0
  for i in 1:nx
    s = s .+ cx[i]*(i+0.0)*x.^(i-1)
  end

  return s
end 
#---------------------------------------------------------------------- 
"""
   function FXY(x,y,c0,cx,cy)

     F(x,y) = c0 + Σcx_i*x^i + Σcy_i*y^i

"""
function FXY(x,y,c0,cx,cy)

  s0 = c0
  sx = FX(x,0.0,cx)
  sy = FX(y,0.0,cy)

  s  = s0 .+ sx .+ sy

  return s
end  
#---------------------------------------------------------------------- 
"""
   function GradFXY(x,y,cx,cy)

     ∇(F(x,y)) = ∇(F(x)) + ∇(F(y))

"""
function GradFXY(x,y,cx,cy)

  sx = GradFX(x,cx)
  sy = GradFX(y,cy)

  s  = sx .+ sy

  return s
end  
#---------------------------------------------------------------------- 

"""
   function FXYZ(x,y,z,c0,cx,cy)

     F(x,y,z) = c0 + Σcx_i*x^i + Σcy_i*y^i + Σcz_i*z^i

"""
function FXYZ(x,y,z,c0,cx,cy,cz)
  nx = length(cx)
  ny = length(cy)
  nz = length(cz)

  s = c0
  for i in 1:nx
    s = s .+ cx[i]*x.^i
  end

  for j in 1:ny
    s = s .+ cy[j]*y.^j
  end  

  for k in 1:nz
    s = s .+ cz[k]*z.^k
  end  
 
  return s
end  
#---------------------------------------------------------------------- 
"""
   function TransFXYZ(x,y,λ,c0,cx,cy)

     x1 = x .- λ*cos(θ)    
     y1 = y .- λ*sin(θ)

     F(x,y) = c0 + Σcx_i*x1^i + Σcy_i*y1^i

"""
function TransFXY(x,y,λ,θ,c0,cx,cy)

  nx = length(cx)
  ny = length(cy)

  x1 = x .- λ*cos(θ)    
  y1 = y .- λ*sin(θ)

  s = c0
  for i in 1:nx
    s = s .+ cx[i]*x1.^i
  end

  for j in 1:ny
    s = s .+ cy[j]*y1.^j
  end  

  return s
end  
#---------------------------------------------------------------------- 
"""
   function RotFXY(x,y,θ,c0,cx,cy)

     Rotate the function with parameters c0,cy,cy about the origin by θ. 

     x1     =  cos(θ)x + sin(θ)y
     y1     = -sin(θ)x + cos(θ)y
     
     F(x,y) = c0 + Σcx_i*x1^i + Σcy_i*y1^i

     Positive θ is an anticlock-wise rotation of the function.
     => Clockwise rotation of the arguments.

"""
function RotFXY(x,y,θ::Float64,c0,cx,cy)

  nx = length(cx)
  ny = length(cy)

  s = c0
  for i in 1:nx
    s = s .+ cx[i]*(cos(θ)*x .+ sin(θ)*y).^i
  end

  for j in 1:ny
    s = s .+ cy[j]*(-sin(θ)*x .+ cos(θ)*y).^j
  end  

  return s
end  
#---------------------------------------------------------------------- 
"""
   function RotLinearFXY(x,y,θ,c0,cx,cy)

     Assuming F(x,y) is a Linear function of both x and y
     Rotate the slope of the function by θ and evaluate the function. 

     F(x,y) = c0 + cy_1*(y + (cx_1/cy_1)*x)

     x1     =  cos(θ)x + sin(θ)y
     y1     = -sin(θ)x + cos(θ)y
     
     F(x,y) = c0 + Σcx_i*x1^i + Σcy_i*y1^i

     Positive θ is an anticlock-wise rotation of the function.
     => Clockwise rotation of the arguments.

"""
function RotLinearFXY(x,y,θ::Float64,c0,cx,cy)

  nx = length(cx)
  ny = length(cy)

  @assert length(cx) == length(cy) == 1 "Non-Linear function provided in RotLinearFXY"

  dm = tan(θ)

  s = c0 .+ cy[1].*(y .+ (cx[1]/cy[1] - dm).*x)

  return s
end  
#---------------------------------------------------------------------- 
"""
   function RotLinearFXY2(x,y,θ,c0,cx,cy)

     Assuming F(x,y) is a Linear function of both x and y
     Rotate the slope of the function by θ and evaluate the function. 

     F(x,y) = c0 + cy_1*(y + (cx_1/cy_1)*x)

     x1     =  cos(θ)x + sin(θ)y
     y1     = -sin(θ)x + cos(θ)y
     
     F(x,y) = c0 + Σcx_i*x1^i + Σcy_i*y1^i

     Positive θ is an anticlock-wise rotation of the function.
     => Clockwise rotation of the arguments.

"""
function RotLinearFXY2(x,y,θ::Float64,c0,cx,cy)

  nx = length(cx)
  ny = length(cy)

  @assert length(cx) == length(cy) == 1 "Non-Linear function provided in RotLinearFXY"

  ϕ  = atan(cx[1],cy[1])
  m  = tan(ϕ - θ)

  s = c0 .+ cy[1].*(y .+ m.*x)

  return s
end  
#---------------------------------------------------------------------- 
"""
   function RotLinearFXY3(x,y,θ,c0,cx,cy)

     Assuming F(x,y) is a Linear function of both x and y
     Rotate the slope of the function by θ and evaluate the function. 

     F(x,y) = c0 + cy_1*(y + (cx_1/cy_1)*x)

     x1     =  cos(θ)x + sin(θ)y
     y1     = -sin(θ)x + cos(θ)y
     
     F(x,y) = c0 + Σcx_i*x1^i + Σcy_i*y1^i

     Positive θ is an anticlock-wise rotation of the function.
     => Clockwise rotation of the arguments.

"""
function RotLinearFXY3(x,y,θ::Float64,c0,cx,cy)

  nx = length(cx)
  ny = length(cy)

  @assert length(cx) == length(cy) == 1 "Non-Linear function provided in RotLinearFXY"

  ϕ  = atan(cx[1],cy[1])
  m  = tan(ϕ - θ)
  R  = sqrt(cy[1]^2 + cx[1]^2)

  s = c0 .+ R*(cos.(ϕ - θ).*y .+ sin.(ϕ - θ).*x)

  return s
end  
#---------------------------------------------------------------------- 


"""
   function RotXYFXYZ(x,y,x0,y0,θ,c0,cx,cy)
      
     x1     = x0 + cos(θ)(x-x0) + sin(θ)(y-y0)
     y1     = y0 - sin(θ)(x-x0) + cos(θ)(y-y0)
     
     F(x,y) = c0 + Σcx_i*x1^i + Σcy_i*y1^i

     Positive θ is an anticlock-wise rotation of the function.
     => Clockwise rotation of the arguments.

"""
function RotXYFXY(x,y,x0::Float64,y0::Float64,θ::Float64,c0,cx,cy)

  nx = length(cx)
  ny = length(cy)

  s = c0
  for i in 1:nx
    s = s .+ cx[i]*(x0 .+ cos(θ)*(x .- x0) .+ sin(θ)*(y .- y0)).^i
  end

  for j in 1:ny
    s = s .+ cy[j]*(y0 .- sin(θ)*(x .- x0) .+ cos(θ)*(y .- y0)).^j
  end  

  return s
end  
#---------------------------------------------------------------------- 
"""
   function RotXYFXYZ(x,y,x0,y0,θ,c0,cx,cy)
      
     x1     = x0 + cos(θ)(x-x0) + sin(θ)(y-y0)
     y1     = y0 - sin(θ)(x-x0) + cos(θ)(y-y0)
     
     F(x,y) = c0 + Σcx_i*x1^i + Σcy_i*y1^i

     Positive θ is an anticlock-wise rotation of the function.
     => Clockwise rotation of the arguments.

"""
function RotXYFXY(x,y,x0::Float64,y0::Float64,θ,c0,cx,cy)

  nx = length(cx)
  ny = length(cy)

  s = c0
  for i in 1:nx
    s = s .+ cx[i]*(x0 .+ cos.(θ).*(x .- x0) .+ sin.(θ).*(y .- y0)).^i
  end

  for j in 1:ny
    s = s .+ cy[j]*(y0 .- sin.(θ).*(x .- x0) .+ cos.(θ).*(y .- y0)).^j
  end  

  return s
end  
#---------------------------------------------------------------------- 

"""
   function TransRotFXYZ(x,y,x0,y0,θ,c0,cx,cy)
      
     x1     = + cos(θ)(x-x0) + sin(θ)(y-y0)
     y1     = - sin(θ)(x-x0) + cos(θ)(y-y0)
     
     F(x,y) = c0 + Σcx_i*x1^i + Σcy_i*y1^i

     Positive θ is an anticlock-wise rotation of the function.
     => Clockwise rotation of the arguments.

"""
function TransRotFXY(x,y,x0::Float64,y0::Float64,θ::Float64,c0,cx,cy)

  nx = length(cx)
  ny = length(cy)

  s = c0
  for i in 1:nx
    s = s .+ cx[i]*(cos(θ)*(x-x0) .+ sin(θ)*(y-y0)).^i
  end

  for j in 1:ny
    s = s .+ cy[j]*(-sin(θ)*(x-x0) .+ cos(θ)*(y-y0)).^j
  end  

  return s
end  
#---------------------------------------------------------------------- 


function TERR_F(x,y,c0,cx,cy)

  n  = length(x)
  nx = length(cx)
  ny = length(cy)

  er = 0.0
  for k in 1:n
    s = c0
    for i in 1:nx
      s = s + cx[i]*x[k]^i
    end

    for j in 1:ny
      s = s + cy[j]*y[k]^j
    end
    er = er + s^2
  end  

  return er
end  
#---------------------------------------------------------------------- 









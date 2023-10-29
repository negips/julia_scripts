#!/bin/julia

#include("NullClines.jl")
#include("ElementOf.jl")

function FXY(x,y,c0,cx,cy)
  nx = length(cx)
  ny = length(cy)

  s = c0
  for i in 1:nx
    s = s + cx[i]*x^i
  end

  for j in 1:ny
    s = s + cy[j]*y^j
  end  

  return s
end  
#---------------------------------------------------------------------- 
function FX(x,c0,cx)
  nx = length(cx)

  s = c0
  for i in 1:nx
    s = s + cx[i]*x^i
  end

  return s
end  
#---------------------------------------------------------------------- 
function FY(y,c0,cy)
  ny = length(cy)

  s = c0
  for i in 1:ny
    s = s + cy[i]*y^i
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

function NullClineFcn(xi,yi,xdx,ydx,xdy,ydy,nx,ny)

  n  = nx + ny + 1      # No of free parameters/Conditions to satisfy
  
  m = length(xi)       # No of incident points
  x = copy(xi)
  y = copy(yi)
  
  # X-Derivative Points
  mx = length(xdx)
  
  # Y-Derivative Points
  my = length(xdy)
  
  # Build the matrix
  N = m + mx + my         # No of constraints to satisfy 
  Mat = zeros(Float64,n,n)
  
  # Incident conditions
  for i in 1:n
    for j in 1:n 
  #   Incident Point contributions      
      for k in 1:m
        if i<=nx+1
          deriv    = (x[k]^(i-1))
        elseif i>nx+1
          deriv    = (y[k]^(i-nx-1))
        end  
       
        if j <= (nx+1)
          Mat[i,j] = Mat[i,j] + (x[k]^(j-1))*deriv
        else
          Mat[i,j] = Mat[i,j] + (y[k]^(j-nx-1))*deriv
        end  
      end     # k
    end       # j
  end         # i
  
  ## Vanishing Gradient contributions df/dx = 0
  Matx = zeros(Float64,n,n)
  for i in 2:nx+1
    for j in 2:nx+1 
  #   x-gradient contributions 
      for k in 1:mx
        deriv     = (i-1)*(xdx[k]^(i-2))
        Matx[i,j] = Matx[i,j] + (j-1)*(xdx[k]^(j-2))*deriv
      end     # k
    end       # j
  end         # i
  
  ## Vanishing Gradient contributions df/dy = 0
  Maty = zeros(Float64,n,n)
  for i in nx+2:n
    for j in nx+2:n 
  #   y-gradient contributions 
      for k in 1:my
        l         = i-(nx+1)
        h         = j-(nx+1)
        deriv     = l*(ydy[k]^(l-1))
        Maty[i,j] = Maty[i,j] + h*(ydy[k]^(h-1))*deriv
      end     # k
    end       # j
  end         # i
  
  M = Mat + Matx + Maty
  
  F = eigen(M)
  j = argmin(F.values)
  r0 = minimum(F.values)
  
  println("Minimum abs Eigenvalue: $r0")
  
  c = F.vectors[:,j]
  
  c0 = c[1]
  cx = c[2:nx+1]
  cy = c[nx+2:nx+ny+1]
  
  return c
end










#!/bin/julia
function element_refine2D(x0::Vector{T},y0::Vector{T},BC0::Vector{String},Par0::AbstractArray{T},vert::Int; eshift=0) where {T<:AbstractFloat}

  @assert length(x0) == 4 "Function only defined for 2D"
  @assert length(y0) == 4 "Function only defined for 2D"

  @assert vert < 5 "Invalid Refinement Vertex"
  @assert vert > 0 "Invalid Refinement Vertex"

  xmid = mean(x0)
  ymid = mean(y0)

  nel  = 3
  ndim = 2
  nc   = 2^ndim
  nf   = 2*ndim

  x   = zeros(T,nc,3)
  y   = zeros(T,nc,3)
  BC  = fill("SYM",nf,nel)
  Par = zeros(T,5,nf,nel) 

  if (vert == 1)    
    # 1st Element
    x[1,1] = x0[1] 
    x[2,1] = 0.5*(x0[1] + x0[2])
    x[3,1] = xmid 
    x[4,1] = 0.5*(x0[1] + x0[4])
    
    y[1,1] = y0[1] 
    y[2,1] = 0.5*(y0[1] + y0[2])
    y[3,1] = ymid 
    y[4,1] = 0.5*(y0[1] + y0[4])

    e       = 1
    BC[1,e] = BC0[1]
    BC[2,e] = "E  "
    BC[3,e] = "E  "
    BC[4,e] = BC0[4]

    f           = 1
    Par[1,f,e]  = Par0[1,1,1]         # Connect to Element No.
    Par[2,f,e]  = Par0[2,1,1]         # Connect on Face No.

    f           = 2
    Par[1,f,e]  = 2 + eshift          # Connect to Element No.
    Par[2,f,e]  = 4                   # Connect on Face No.

    f           = 3
    Par[1,f,e]  = 3 + eshift                   # Connect to Element No.
    Par[2,f,e]  = 3                   # Connect on Face No.

    f           = 4
    Par[1,f,e]  = Par0[1,4,1]         # Connect to Element No.
    Par[2,f,e]  = Par0[2,4,1]         # Connect on Face No.
   
    # Second Element
    x[1,2] = 0.5*(x0[1] + x0[2]) 
    x[2,2] = x0[2] 
    x[3,2] = x0[3]
    x[4,2] = xmid
    
    y[1,2] = 0.5*(y0[1] + y0[2]) 
    y[2,2] = y0[2] 
    y[3,2] = y0[3] 
    y[4,2] = ymid

    e       = 2
    BC[1,e] = BC0[1]
    BC[2,e] = BC0[2]
    BC[3,e] = "E  "
    BC[4,e] = "E  "

    f           = 1
    Par[1,f,e]  = Par0[1,1,1]         # Connect to Element No.
    Par[2,f,e]  = Par0[2,1,1]         # Connect on Face No.

    f           = 2
    Par[1,f,e]  = Par0[1,2,1]         # Connect to Element No.
    Par[2,f,e]  = Par0[2,2,1]         # Connect on Face No.

    f           = 3
    Par[1,f,e]  = 3 + eshift                   # Connect to Element No.
    Par[2,f,e]  = 4                   # Connect on Face No.

    f           = 4
    Par[1,f,e]  = 1 + eshift                   # Connect to Element No.
    Par[2,f,e]  = 2                   # Connect on Face No.


    # Third Element
    x[1,3] = x0[3] 
    x[2,3] = x0[4]
    x[3,3] = 0.5*(x0[4] + x0[1])
    x[4,3] = xmid
    
    y[1,3] = y0[3] 
    y[2,3] = y0[4]
    y[3,3] = 0.5*(y0[4] + y0[1])
    y[4,3] = ymid

    e       = 3
    BC[1,e] = BC0[3]
    BC[2,e] = BC0[4]
    BC[3,e] = "E  "
    BC[4,e] = "E  "

    f           = 1
    Par[1,f,e]  = Par0[1,3,1]         # Connect to Element No.
    Par[2,f,e]  = Par0[2,3,1]         # Connect on Face No.

    f           = 2
    Par[1,f,e]  = Par0[1,2,1]         # Connect to Element No.
    Par[2,f,e]  = Par0[2,2,1]         # Connect on Face No.

    f           = 3
    Par[1,f,e]  = 3 + eshift                   # Connect to Element No.
    Par[2,f,e]  = 4                   # Connect on Face No.

    f           = 4
    Par[1,f,e]  = 1                   # Connect to Element No.
    Par[2,f,e]  = 2                   # Connect on Face No.


  elseif (vert == 2)
    # 1st Element
    x[1,1] = x0[1] 
    x[2,1] = 0.5*(x0[1] + x0[2])
    x[3,1] = xmid 
    x[4,1] = x0[4]
    
    y[1,1] = y0[1] 
    y[2,1] = 0.5*(y0[1] + y0[2])
    y[3,1] = ymid 
    y[4,1] = y0[4]

    e       = 1
    BC[1,e] = BC0[1]
    BC[2,e] = "E  "
    BC[3,e] = "E  "
    BC[4,e] = BC0[4]

    f           = 1
    Par[1,f,e]  = Par0[1,1,1]         # Connect to Element No.
    Par[2,f,e]  = Par0[2,1,1]         # Connect on Face No.

    f           = 2
    Par[1,f,e]  = 2 + eshift                   # Connect to Element No.
    Par[2,f,e]  = 4                   # Connect on Face No.

    f           = 3
    Par[1,f,e]  = 3 + eshift                   # Connect to Element No.
    Par[2,f,e]  = 3                   # Connect on Face No.

    f           = 4
    Par[1,f,e]  = Par0[1,4,1]         # Connect to Element No.
    Par[2,f,e]  = Par0[2,4,1]         # Connect on Face No.

    # Second Element
    x[1,2] = 0.5*(x0[1] + x0[2]) 
    x[2,2] = x0[2] 
    x[3,2] = 0.5*(x0[2] + x0[3])
    x[4,2] = xmid
    
    y[1,2] = 0.5*(y0[1] + y0[2]) 
    y[2,2] = y0[2] 
    y[3,2] = 0.5*(y0[2] + y0[3]) 
    y[4,2] = ymid

    e       = 2
    BC[1,e] = BC0[1]
    BC[2,e] = BC0[2]
    BC[3,e] = "E  "
    BC[4,e] = "E  "

    f           = 1
    Par[1,f,e]  = Par0[1,1,1]         # Connect to Element No.
    Par[2,f,e]  = Par0[2,1,1]         # Connect on Face No.

    f           = 2
    Par[1,f,e]  = Par0[1,2,1]         # Connect to Element No.
    Par[2,f,e]  = Par0[2,2,1]         # Connect on Face No.

    f           = 3
    Par[1,f,e]  = 3 + eshift                   # Connect to Element No.
    Par[2,f,e]  = 4                   # Connect on Face No.

    f           = 4
    Par[1,f,e]  = 1 + eshift                   # Connect to Element No.
    Par[2,f,e]  = 2                   # Connect on Face No.

    # Third Element
    x[1,3] = 0.5*(x0[2] + x0[3]) 
    x[2,3] = x0[3]
    x[3,3] = x0[4] 
    x[4,3] = xmid 
    
    y[1,3] = 0.5*(y0[2] + y0[3]) 
    y[2,3] = y0[3]  
    y[3,3] = y0[4]
    y[4,3] = ymid 

    e       = 3
    BC[1,e] = BC0[2]
    BC[2,e] = BC0[3]
    BC[3,e] = "E  "
    BC[4,e] = "E  "

    f           = 1
    Par[1,f,e]  = Par0[1,2,1]         # Connect to Element No.
    Par[2,f,e]  = Par0[2,2,1]         # Connect on Face No.

    f           = 2
    Par[1,f,e]  = Par0[1,3,1]         # Connect to Element No.
    Par[2,f,e]  = Par0[2,3,1]         # Connect on Face No.

    f           = 3
    Par[1,f,e]  = 1 + eshift                   # Connect to Element No.
    Par[2,f,e]  = 3                   # Connect on Face No.

    f           = 4
    Par[1,f,e]  = 2 + eshift                   # Connect to Element No.
    Par[2,f,e]  = 3                   # Connect on Face No.

  elseif (vert == 3)
    # 1st Element
    x[1,1] = x0[1] 
    x[2,1] = x0[2]
    x[3,1] = 0.5*(x0[2] + x0[3])
    x[4,1] = xmid
    
    y[1,1] = y0[1] 
    y[2,1] = y0[2]
    y[3,1] = 0.5*(y0[2] + y0[3])
    y[4,1] = ymid

    e       = 1
    BC[1,e] = BC0[1]
    BC[2,e] = "E  "
    BC[3,e] = "E  "
    BC[4,e] = BC0[4]

    f           = 1
    Par[1,f,e]  = Par0[1,1,1]         # Connect to Element No.
    Par[2,f,e]  = Par0[2,1,1]         # Connect on Face No.

    f           = 2
    Par[1,f,e]  = Par0[1,2,1]         # Connect to Element No.
    Par[2,f,e]  = Par0[2,2,1]         # Connect on Face No.

    f           = 3
    Par[1,f,e]  = 2 + eshift                   # Connect to Element No.
    Par[2,f,e]  = 3                   # Connect on Face No.

    f           = 4
    Par[1,f,e]  = 3 + eshift                   # Connect to Element No.
    Par[2,f,e]  = 3                   # Connect on Face No.

    # Second Element
    x[1,2] = 0.5*(x0[2] + x0[3]) 
    x[2,2] = x0[3] 
    x[3,2] = 0.5*(x0[3] + x0[4])
    x[4,2] = xmid
    
    y[1,2] = 0.5*(y0[2] + y0[3]) 
    y[2,2] = y0[3] 
    y[3,2] = 0.5*(y0[3] + y0[4]) 
    y[4,2] = ymid

    e       = 2
    BC[1,e] = BC0[2]
    BC[2,e] = BC0[3]
    BC[3,e] = "E  "
    BC[4,e] = "E  "

    f           = 1
    Par[1,f,e]  = Par0[1,2,1]         # Connect to Element No.
    Par[2,f,e]  = Par0[2,2,1]         # Connect on Face No.

    f           = 2
    Par[1,f,e]  = Par0[1,3,1]         # Connect to Element No.
    Par[2,f,e]  = Par0[2,3,1]         # Connect on Face No.

    f           = 3
    Par[1,f,e]  = 3 + eshift                   # Connect to Element No.
    Par[2,f,e]  = 4                   # Connect on Face No.

    f           = 4
    Par[1,f,e]  = 1 + eshift                   # Connect to Element No.
    Par[2,f,e]  = 3                   # Connect on Face No.

    # Third Element
    x[1,3] = 0.5*(x0[3] + x0[4]) 
    x[2,3] = x0[4]
    x[3,3] = x0[1] 
    x[4,3] = xmid 
    
    y[1,3] = 0.5*(y0[3] + y0[4]) 
    y[2,3] = y0[4]  
    y[3,3] = y0[1]
    y[4,3] = ymid 

    e       = 3
    BC[1,e] = BC0[3]
    BC[2,e] = BC0[4]
    BC[3,e] = "E  "
    BC[4,e] = "E  "

    f           = 1
    Par[1,f,e]  = Par0[1,3,1]         # Connect to Element No.
    Par[2,f,e]  = Par0[2,3,1]         # Connect on Face No.

    f           = 2
    Par[1,f,e]  = Par0[1,4,1]         # Connect to Element No.
    Par[2,f,e]  = Par0[2,4,1]         # Connect on Face No.

    f           = 3
    Par[1,f,e]  = 1 + eshift                   # Connect to Element No.
    Par[2,f,e]  = 4                   # Connect on Face No.

    f           = 4
    Par[1,f,e]  = 2 + eshift                   # Connect to Element No.
    Par[2,f,e]  = 3                   # Connect on Face No.

  elseif (vert == 4)
    # 1st Element
    x[1,1] = x0[1] 
    x[2,1] = x0[2]
    x[3,1] = xmid
    x[4,1] = 0.5*(x0[1] + x0[4])
    
    y[1,1] = y0[1] 
    y[2,1] = y0[2]
    y[3,1] = ymid
    y[4,1] = 0.5*(y0[1] + y0[4])

    e       = 1
    BC[1,e] = BC0[1]
    BC[2,e] = "E  "
    BC[3,e] = "E  "
    BC[4,e] = BC0[4]

    f           = 1
    Par[1,f,e]  = Par0[1,1,1]         # Connect to Element No.
    Par[2,f,e]  = Par0[2,1,1]         # Connect on Face No.

    f           = 2
    Par[1,f,e]  = 2 + eshift                   # Connect to Element No.
    Par[2,f,e]  = 4                   # Connect on Face No.

    f           = 3
    Par[1,f,e]  = 2 + eshift                   # Connect to Element No.
    Par[2,f,e]  = 3                   # Connect on Face No.

    f           = 4
    Par[1,f,e]  = Par0[1,4,1]         # Connect to Element No.
    Par[2,f,e]  = Par0[2,4,1]         # Connect on Face No.

    # Second Element
    x[1,2] = x0[2] 
    x[2,2] = x0[3] 
    x[3,2] = 0.5*(x0[3] + x0[4])
    x[4,2] = xmid
    
    y[1,2] = y0[2]
    y[2,2] = y0[3] 
    y[3,2] = 0.5*(y0[3] + y0[4]) 
    y[4,2] = ymid

    e       = 2
    BC[1,e] = BC0[2]
    BC[2,e] = BC0[3]
    BC[3,e] = "E  "
    BC[4,e] = "E  "

    f           = 1
    Par[1,f,e]  = Par0[1,2,1]         # Connect to Element No.
    Par[2,f,e]  = Par0[2,2,1]         # Connect on Face No.

    f           = 2
    Par[1,f,e]  = Par0[1,3,1]         # Connect to Element No.
    Par[2,f,e]  = Par0[2,3,1]         # Connect on Face No.

    f           = 3
    Par[1,f,e]  = 3 + eshift                   # Connect to Element No.
    Par[2,f,e]  = 4                   # Connect on Face No.

    f           = 4
    Par[1,f,e]  = 1 + eshift                   # Connect to Element No.
    Par[2,f,e]  = 2                   # Connect on Face No.

    # Third Element
    x[1,3] = 0.5*(x0[3] + x0[4]) 
    x[2,3] = x0[4]
    x[3,3] = 0.5*(x0[4] + x0[1]) 
    x[4,3] = xmid 
    
    y[1,3] = 0.5*(y0[3] + y0[4]) 
    y[2,3] = y0[4]  
    y[3,3] = 0.5*(y0[4] + y0[1])
    y[4,3] = ymid 

    e       = 3
    BC[1,e] = BC0[3]
    BC[2,e] = BC0[4]
    BC[3,e] = "E  "
    BC[4,e] = "E  "

    f           = 1
    Par[1,f,e]  = Par0[1,3,1]         # Connect to Element No.
    Par[2,f,e]  = Par0[2,3,1]         # Connect on Face No.

    f           = 2
    Par[1,f,e]  = Par0[1,4,1]         # Connect to Element No.
    Par[2,f,e]  = Par0[2,4,1]         # Connect on Face No.

    f           = 3
    Par[1,f,e]  = 1 + eshift                   # Connect to Element No.
    Par[2,f,e]  = 3                   # Connect on Face No.

    f           = 4
    Par[1,f,e]  = 2 + eshift                   # Connect to Element No.
    Par[2,f,e]  = 3                   # Connect on Face No.

  end  

  return x,y,BC,Par
end
#---------------------------------------------------------------------- 
function twoelements_refine2D(x0::Matrix{T},y0::Matrix{T},BC0::Matrix{String},Par0::Array{T},verts::Vector{Int64}) where {T<:AbstractFloat}

  nel2  = 6
  ndim  = 2
  nf    = 2*ndim
  nc    = 2^ndim

  x     = zeros(T,nc,nel2)
  y     = zeros(T,nc,nel2)
  BC    = fill("E  ",nf,nel2)
  Par   = zeros(T,5,nf,nel2)
  es1   = 0
  es2   = 3
  
  x[:,1:3],y[:,1:3],BC[:,1:3],Par[:,:,1:3] = element_refine2D(x0[:,1],y0[:,1],BC0[:,1],Par0[:,:,1],verts[1];eshift=es1)
  x[:,4:6],y[:,4:6],BC[:,4:6],Par[:,:,4:6] = element_refine2D(x0[:,2],y0[:,2],BC0[:,2],Par0[:,:,2],verts[2];eshift=es2)

  if (verts[1] == 2 && verts[2] == 1)
    e = 2
    f = 2
    Par[1,f,e] = es2+1  # Element No. 
    Par[2,f,e] = 4      # On Face

    e = 3
    f = 4
    Par[1,f,e] = es2+3  # Element No. 
    Par[2,f,e] = 2      # On Face

    e = 4
    f = 4
    Par[1,f,e] = es1+2  # Element No. 
    Par[2,f,e] = 2      # On Face

    e = 6
    f = 2
    Par[1,f,e] = es1+3  # Element No. 
    Par[2,f,e] = 4      # On Face

  elseif (verts[1] == 1 && verts[2] == 2) 
    e = 1
    f = 4
    Par[1,f,e] = es2+2  # Element No. 
    Par[2,f,e] = 2      # On Face

    e = 3
    f = 2
    Par[1,f,e] = es2+3  # Element No. 
    Par[2,f,e] = 4      # On Face

    e = 5
    f = 2
    Par[1,f,e] = es1+1  # Element No. 
    Par[2,f,e] = 4      # On Face

    e = 6
    f = 4
    Par[1,f,e] = es1+3  # Element No. 
    Par[2,f,e] = 2      # On Face

  elseif (verts[1] == 3 && verts[2] == 4)
    e = 1
    f = 2
    Par[1,f,e] = es2+1  # Element No. 
    Par[2,f,e] = 4      # On Face

    e = 2
    f = 1
    Par[1,f,e] = es2+3  # Element No. 
    Par[2,f,e] = 2      # On Face

    e = 4
    f = 4
    Par[1,f,e] = es1+1  # Element No. 
    Par[2,f,e] = 2      # On Face

    e = 6
    f = 2
    Par[1,f,e] = es1+2  # Element No. 
    Par[2,f,e] = 1      # On Face

  elseif (verts[1] == 4 && verts[2] == 3)
    e = 1
    f = 4
    Par[1,f,e] = es2+1  # Element No. 
    Par[2,f,e] = 2      # On Face

    e = 3
    f = 2
    Par[1,f,e] = es2+2  # Element No. 
    Par[2,f,e] = 1      # On Face

    e = 4
    f = 2
    Par[1,f,e] = es1+1  # Element No. 
    Par[2,f,e] = 4      # On Face

    e = 5
    f = 1
    Par[1,f,e] = es1+3  # Element No. 
    Par[2,f,e] = 2      # On Face

  end  


  return x,y,BC,Par
end
#---------------------------------------------------------------------- 
function edge_refine2D(x0::Matrix{T},y0::Matrix{T},BC0::Matrix{String},Par0::Array{T},verts::Vector{Int64}) where {T<:AbstractFloat}

  nel2  = 9
  ndim  = 2
  nf    = 2*ndim
  nc    = 2^ndim

  x     = zeros(T,nc,nel2)
  y     = zeros(T,nc,nel2)
  BC    = fill("E  ",nf,nel2)
  Par   = zeros(T,2,nf,nel2)
  
  verts = [2; 1; 4]
 
  x[:,1:6],y[:,1:6],BC[:,1:6],Par[:,:,1:6] = twoelements_refine2D(x0,y0,BC0,Par0,verts)

  es3 = 6
  x[:,7:9],y[:,7:9],BC[:,7:9],Par[:,:,7:9] = element_refine2D(x0[:,3],y0[:,3],BC0[:,3],Par0[:,:,3],verts[3];eshift=es3)


  return x,y,BC,Par
end
#---------------------------------------------------------------------- 












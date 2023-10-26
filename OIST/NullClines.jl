function NullClines(fxy,xi,yr0,yr1,nsteps,dτ)

  x    = xi
  y    = xi
  
  fx(z) = fxy(z,y)
  fy(z) = fxy(x,z)
  
  # Initial point
  x0 = xi
  x  = x0
  yr = find_zeros(fy,yr0,yr1,verbose=true)
  
  nroots   = length(yr)
  println("Nroots = $nroots")

  if nroots == 0

    display("No roots found within range ($yr0,$yr1) for xi=$xi")
    return xi,xi

  end  

  y  = yr[1]
  froots_y = zeros(nsteps,nroots)
  froots_x = zeros(nsteps,nroots)

  stepmax = abs(dτ)
  
  # Find null-cline
  for n in 1:nroots
  
    x   = x0
    yn  = find_zeros(fy,yr0,yr1);
    x   = x0 + dτ
    yn1 = find_zeros(fy,yr0,yr1);

    dy = yn1[n] - yn[n]
    y  = yn[n]
    x  = x0
    dx = dτ
  
    for i in 1:nsteps
    
      if abs(dy) > abs(dx)
        y  = y + dy
        x1 = find_zero(fx,x,atol=1.0e-8,verbose=false);
        dx = x1 - x
        x  = x1
      else
        x  = x + dx
        y1 = find_zero(fy,y,atol=1.0e-08,verbose=false);
        dy = y1 - y
        y  = y1
      end
    
      froots_y[i,n] = y
      froots_x[i,n] = x
    end
  end  

  return froots_x,froots_y
end
#---------------------------------------------------------------------- 


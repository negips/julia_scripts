#!/bin/julia

#println("Calculate the null clines for Model Equations in H. Meinhardt (1995) The algorithmic beauty of sea shells")
function Meinhardt_Nullclines(EQN,Mein)

#using Roots

#include("Meinhardt.jl")

#EQN = "2.4"

if (EQN == "2.1")
#  Mein(x,y,s)   = Meinhardt_21(x,y,s)
  arange0_a  = -40.0
  arange1_a  = +40.0
  b0_a       = 25.0
  τa         = -1.0e-3

  arange0_b  = -40.0
  arange1_b  = +40.0
  b0_b       = 25.0
  τb         = -1.0e-3
elseif (EQN == "2.4")
#  Mein(x,y,s)   = Meinhardt_24(x,y,s)
  arange0_a  = -80.0
  arange1_a  = +80.0
  b0_a       =  15.0
  τa         = -5.0e-4

  arange0_b  = -80.0
  arange1_b  = +80.0
  b0_b       =  0.45
  τb         = 1.0e-4
 
elseif (EQN == "2.5")
#  Mein(x,y,s)   = Meinhardt_25(x,y,s)
  arange0_a  = -40.0
  arange1_a  = +40.0
  b0_a       = 25.0
  τa         = -1.0e-3

  arange0_b  = -40.0
  arange1_b  = +40.0
  b0_b       = -25.0
  τb         = 1.0e-3
elseif (EQN == "3.1")
#  Mein(x,y,s)   = Meinhardt_31(x,y,s)
  arange0_a  = -40.0
  arange1_a  = +40.0
  b0_a       = 25.0
  τa         = -1.0e-3

  arange0_b  = -40.0
  arange1_b  = +40.0
  b0_b       = +25.0
  τb         = -1.0e-3
elseif (EQN == "6.1")
#  Mein(x,y,s)   = Meinhardt_31(x,y,s)
  arange0_a  = -50.0
  arange1_a  = +50.0
  b0_a       = 10.0
  τa         = -1.0e-3

  arange0_b  = -50.0
  arange1_b  = +50.0
  b0_b       = +10.0
  τb         = -1.0e-3
elseif (EQN == "1987_1")
#  Mein(x,y,s)   = Meinhardt_31(x,y,s)
  arange0_a  = -70.0
  arange1_a  = +70.0
  b0_a       =  5.0
  τa         = -1.0e-4

  arange0_b  = -10.0
  arange1_b  = +10.0
  b0_b       = +5.0
  τb         = -1.0e-3
elseif (EQN == "1987_2")
#  Mein(x,y,s)   = Meinhardt_31(x,y,s)
  arange0_a  = -10.0
  arange1_a  = +10.0
  b0_a       =  35.0
  τa         = -1.0e-3

  arange0_b  = -10.0
  arange1_b  = +10.0
  b0_b       = +10.0
  τb         = -1.0e-3
elseif (EQN == "1987_2_branching")
#  Mein(x,y,s)   = Meinhardt_31(x,y,s)
  arange0_a  = -10.0
  arange1_a  = +10.0
  b0_a       = +15.0
  τa         = -1.0e-3

  arange0_b  = -10.0
  arange1_b  = +10.0
  b0_b       =  1.0
  τb         = -4.0e-4
 
else
  display("$EQN not defined.")
end

s  = 1

Adotxy(x,y) = Mein(x,y)[1]
Adotx(x) = Adotxy(x,b)
Adoty(y) = Adotxy(a,y)

nsteps = 50000
dτ     = τa
aroots_a,aroots_b = find_nullcline(Adotxy,b0_a,b0_a,arange0_a,arange1_a,nsteps,dτ)

Bdotxy(x,y) = Mein(x,y)[2]
Bdotx(x) = Bdotxy(x,b)
Bdoty(y) = Bdotxy(a,y)
dτ     = τb
broots_a,broots_b = find_nullcline(Bdotxy,b0_b,b0_b,arange0_b,arange1_b,nsteps,dτ)

return aroots_a,aroots_b,broots_a,broots_b
end


#---------------------------------------------------------------------- 

function find_nullcline(fxy,ai,bi,ar0,ar1,nsteps,dτ)

  b    = bi
  a    = ai
  
  fx(x) = fxy(x,b)
  fy(y) = fxy(a,y)
  
  # Initial point
  b0 = bi
  b  = b0
  ar = find_zeros(fx,ar0,ar1)
  
  nroots   = length(ar)
  println("Nroots = $nroots")

  if nroots == 0

    display("No roots found within range ($ar0,$ar1)")
    return ai,bi

  end  

  a  = ar[1]
  aroots_a = zeros(nsteps,nroots)
  aroots_b = zeros(nsteps,nroots)

  stepmax = abs(dτ)
  
  # Find Activator null-cline
  for n in 1:nroots
  
    b   = b0
    an  = find_zeros(fx,ar0,ar1);
    b   = b0 + dτ
    an1 = find_zeros(fx,ar0,ar1);

    da = an1[n] - an[n]
    a  = an[n]
    b  = b0
    db = dτ
  
    for i in 1:nsteps
    
      if abs(da) > abs(db)
#        if abs(da) > stepmax
#          da = sign(da)*dτ
#        end  
        a  = a + da
        b1 = find_zero(fy,b,atol=1.0e-8,verbose=false);
        db = b1 - b
        b  = b1
      else
#        if abs(db) > stepmax
#          db = sign(db)*dτ
#        end  
        b  = b + db
        a1 = find_zero(fx,a,atol=1.0e-08,verbose=false);
        da = a1 - a
        a  = a1
      end
    
      aroots_a[i,n] = a
      aroots_b[i,n] = b
    end
  end  

  return aroots_a,aroots_b
end
#---------------------------------------------------------------------- 








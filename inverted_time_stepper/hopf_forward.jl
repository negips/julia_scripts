println("Time integrate the standard Hopf normal form")

using PyPlot
using LinearAlgebra

include("GetBDF.jl")
include("GetEXT.jl")

close("all")

μ   = 1.0

x  = 0.0;
y  = 0.0

# Equations:
#
# ̇x = μx + y
# ̇y = -x + μy -x²y
#

x0 = 0.01
y0 = 0.0

bdford = 2
extord = 2


nstep = 100000
dt    = 0.0001

xhist = zeros(Float64,nstep+1)
yhist = zeros(Float64,nstep+1)


x  = x0;
y  = y0

xlag = zeros(Float64,2)
ylag = zeros(Float64,2)

fxlag = zeros(Float64,2)
fylag = zeros(Float64,2)

xhist[1] = x
yhist[1] = y

for i in 1:nstep
   global x,y,xlag,ylag,fxlag,fylag

   bdfk = zeros(Float64,4)
   extk = zeros(Float64,3)

   if i==1
     bord = 1
     eord = 1
   elseif i==2  
     bord = 2
     eord = 2
   else
     bord = bdford
     eord = extord
   end

   if bord>bdford
     bord = bdford
   end
   if eord>extord
     eord = extord
   end  

   GetBDF!(bdfk,bord)
   GetEXT!(extk,eord)

   rhsx = -(bdfk[2]*x + bdfk[3]*xlag[1] + bdfk[4]*xlag[2])/dt
   rhsy = -(bdfk[2]*y + bdfk[3]*ylag[1] + bdfk[4]*ylag[2])/dt

   xlag[1] = x
   xlag[2] = xlag[1]

   ylag[1] = y
   ylag[2] = ylag[1]

   fx   = 0
   fy   = -x*x*y

   rhsx = rhsx + extk[1]*fx + extk[2]*fxlag[1] + extk[3]*fxlag[2]
   rhsy = rhsy + extk[1]*fy + extk[2]*fylag[1] + extk[3]*fylag[2]

   fxlag[1] = fx
   fxlag[2] = fxlag[1]

   fylag[1] = fy
   fylag[2] = fylag[1]

   H    = [bdfk[1]/dt-μ -1;
           1            bdfk[1]/dt-μ]

   rhs  = [rhsx; rhsy]

   s    = inv(H)*rhs
   x    = s[1]
   y    = s[2]

#   x    = rhsx/(bdfk[1]/dt - μ)
#   y    = rhsy/(bdfk[1]/dt - μ)

#   println(x)
#   println(y)

   xhist[i+1] = x
   yhist[i+1] = y

end


plot(xhist,yhist)










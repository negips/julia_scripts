module mB

   include("mA.jl")

   using .mA: a as b

   global B::Array

   export fillb!

   function fillb!(r::Float64)

       fill!(b,r)

       return nothing
   end

   function create_B(r::Float64,tup::Tuple)
       
       global B = fill(r,tup)

       return nothing
   end


end


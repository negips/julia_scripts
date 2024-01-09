module mA

   const n = 10
   a = Vector{Float64}(undef,n)
  
   # A::Array

   export filla!, create_A

   function filla!(r::Float64)

       fill!(a,r)

       return nothing
   end

   function create_A(r::Float64,tup::Tuple)
       
       global A = fill(r,tup)

       return nothing
   end


end


#     Add Definition of Constructors here
#---------------------------------------------------------------------- 
"""
      function SEM_Input(N::Int,Nd::Int,nel::Int,xs::Float64,xe::Float64)

      N     - Polynomial Order
      Nd    - Dealiasing Order
      nel   - No. of elements
      xs    - Start of domain
      xe    - End of domain
      Dtype - Datatype (Float64,ComplexF64,...)

"""
function SEM_Input(N::Int,Nd::Int,nel::Int,xs::Float64,xe::Float64;Dtype=Float64)

  lx1       = N+1
  lxd       = Nd+1
  nnodes    = nel+1
  xc        = Vector(range(xs,stop=xe,length=nnodes))

  return SEM_Input(N,lx1,Nd,lxd,nel,xs,xe,xc,Dtype)
end
#----------------------------------------------------------------------       

 


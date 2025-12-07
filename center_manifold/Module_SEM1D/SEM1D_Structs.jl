#     Add Definition of new structures here
#---------------------------------------------------------------------- 
"""
      struct SEMInput

      Has fields:

      N     - Polynomial Order
      lx1   - No of Points.
      Nd    - Dealiasing Order
      lxd   - No. of Dealiasing points
      nel   - No. of elements
      xs    - Start of domain
      xe    - End of domain
      xc    - Vector of nodes
      lbc   - left boundary condition: true => Dirichlet
      rbc   - right boundary condition true => Dirichlet
      Dtype - Datatype

"""
struct SEMInput <: AbstractParams

      N::Int
      lx1::Int
      Nd::Int
      lxd::Int
      nel::Int
      xs::Float64
      xe::Float64
      xc::Vector{Float64}
      lbc::Bool
      rbc::Bool
      Dtype::DataType

end
#----------------------------------------------------------------------       

"""
      struct SEMGeoMat{T}

      Has fields:

      xm1         - Jacobi nodes (GLL/GL)
      xrm1        - dx/dr
      rxm1        - dr/dx
      jacm1       - Jacobian Matrix
      jacmi       - Jacobian Inverse
      bm1         - Mass (Vector)
      Gradx       - Gradient
      Intpm1d     - Interpolation lx1 -> lxd
      Gradxd      - Gradient on lxd
      bm1d        - Mass (Vector) on lxd 
      Bintpd      - Integral on lxd (bm1d*intpm1d)' 
      Convd       - Convection operator (with dealiasing) 
      WLap        - Weak Laplacian
      Lap         - Laplacian w/o integration by parts

"""
struct SEMGeoMat{T<:Number} <: AbstractGeometry

      xm1;
      xrm1;
      rxm1;
      jacm1;
      jacmi;
      bm1;
      Gradx;
      Intpm1d;
      Gradxd;
      bm1d;
      Bintpd;
      Convd;
      WLap;
      Lap;
end
#----------------------------------------------------------------------





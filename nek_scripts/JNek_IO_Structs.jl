#     Add Definition of new structures here
"""
      mutable struct TwoTensorField{T}

        Has fields:

        dims::Int             # Dimension of the Tensor
        p::Vector{Int}        # Polynomial orders in the different directions
        nel::Int              # No of Tensor elements
        field::Array{T}       # Tensor Field

"""
      mutable struct TwoTensorField{T <: Number} <: AbstractTensorField{T}
        dims::Int
        p::Vector{Int}
        nel::Int
        field::Array{T}
      end
#---------------------------------------------------------------------- 
#     Corresponding Function for Initialization
"""
        function TwoTensorField(field::Array)
          
        Constructor for the TwoTensorField struct. 

"""
     function TwoTensorField(field::Array{T}) where T <: Number

        el   = eltype(field)
        s    = collect(Int,size(field))
        @assert length(s) == 3
        
        dims = 2
        p    = s[1:2]
        nel  = s[3]

        return TwoTensorField{el}(dims,p,nel,field)
      end        

#---------------------------------------------------------------------- 
      mutable struct NekField{T} #where T <: AbstractFloat
#       4/8 Byte Fields        
        hdr::String
        version::String
        wdsize::Int
        lx1::Int
        ly1::Int
        lz1::Int
        nel::Int
        nelgt::Int
        time::T
        istep::Int
        fid0::Int
        nfileo::Int
        rdcode::String
        p0th::T
        ifprmesh::Bool
        glnum::Vector{Int}
        x::Array{T}
        y::Array{T}
        z::Array{T}
        u::Array{T}
        v::Array{T}
        w::Array{T}
        p::Array{T}
        t::Array{T}

      end
#---------------------------------------------------------------------- 

      mutable struct Re2Field{T} #where T <: AbstractFloat

#       8 Byte Fields
        wdsize::Int
        hdr::String
        version::String
        nelgt::Int
        ldimr::Int
        nelgv::Int
        xc::Array{T}
        yc::Array{T}
        zc::Array{T}
        ncurve::Int
        curveieg::Vector{Int}
        curveiside::Vector{Int}
        curveparam::Array{T}
        curvetype::Vector{String}
        cbl::Array{String}
        bl::Array{T}

      end       
#---------------------------------------------------------------------- 
     mutable struct ma2Hdr

#       .ma2 data
#        hdr::String
        version::String
        nel::Int 
        nactive::Int 
        depth::Int 
        d2::Int 
        npts::Int 
        nrank::Int 
        noutflow::Int 

      end       

#----------------------------------------------------------------------  
     mutable struct ma2Field

#       .ma2 data
#        hdr::String
        pmap::Vector{Int}
        vmap::Array{Int}

      end       

#----------------------------------------------------------------------  

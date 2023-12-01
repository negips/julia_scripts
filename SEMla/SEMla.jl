#     Port for reader_par.f
#     Author:     Prabal Negi
#

      module SEMla

      export read_re2_hdr

      mutable struct NekField8
#       8 Byte Fields        
        hdr::String
        version::String
        wdsize::Int
        lx1::Int
        ly1::Int
        lz1::Int
        nel::Int
        nelgt::Int
        time::Float32
        istep::Int
        fid0::Int
        nfileo::Int
        rdcode::String
        p0th::Float32
        ifprmesh::Bool
        glnum::Vector{Int}
        x::Array{Float32}
        y::Array{Float32}
        z::Array{Float32}
        u::Array{Float32}
        v::Array{Float32}
        w::Array{Float32}
        p::Array{Float32}
        t::Array{Float32}

      end
#---------------------------------------------------------------------- 
      mutable struct NekField16
#       16 Byte fields        
        hdr::String
        version::String
        wdsize::Int
        lx1::Int
        ly1::Int
        lz1::Int
        nel::Int
        nelgt::Int
        time::Float64
        istep::Int
        fid0::Int
        nfileo::Int
        rdcode::String
        p0th::Float64
        ifprmesh::Bool
        glnum::Vector{Int}
        x::Array{Float64}
        y::Array{Float64}
        z::Array{Float64}
        u::Array{Float64}
        v::Array{Float64}
        w::Array{Float64}
        p::Array{Float64}
        t::Array{Float64}

      end  
#---------------------------------------------------------------------- 
     mutable struct Re2Field8

#       8 Byte Fields
        wdsize::Int
        hdr::String
        version::String
        nelgt::Int
        ldimr::Int
        nelgv::Int
        xc::Array{Float32}
        yc::Array{Float32}
        zc::Array{Float32}
        ncurve::Int
        curveieg::Vector{Int}
        curveiside::Vector{Int}
        curveparam::Array{Float32}
        curvetype::Vector{String}
        cbl::Array{String}
        bl::Array{Float32}

      end       
#---------------------------------------------------------------------- 
     mutable struct Re2Field16

#       16 Byte Fields
        wdsize::Int
        hdr::String
        version::String
        nelgt::Int
        ldimr::Int
        nelgv::Int
        xc::Array{Float64}
        yc::Array{Float64}
        zc::Array{Float64}
        ncurve::Int
        curveieg::Vector{Int}
        curveiside::Vector{Int}
        curveparam::Array{Float64}
        curvetype::Vector{String}
        cbl::Array{String}
        bl::Array{Float64}

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

      function __init__()


        return nothing
      end 

#----------------------------------------------------------------------

      function init()

#       Can't reinitialize MPI        

        return nothing
      end 

#----------------------------------------------------------------------

#---------------------------------------------------------------------- 
      end   # Module JNek_IO_MPI











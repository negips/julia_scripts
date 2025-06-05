#     Add Definition of new structures here
#---------------------------------------------------------------------- 
      mutable struct NekFldHdr{T} 

#       4/8 Byte Fields        
#        hdr::String
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

      end
#---------------------------------------------------------------------- 
      mutable struct NekField{T} #where T <: AbstractFloat
#       4/8 Byte Fields        
        hdr::NekFldHdr
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
#      mutable struct NekField{T} #where T <: AbstractFloat
##       4/8 Byte Fields        
#        hdr::String
#        version::String
#        wdsize::Int
#        lx1::Int
#        ly1::Int
#        lz1::Int
#        nel::Int
#        nelgt::Int
#        time::T
#        istep::Int
#        fid0::Int
#        nfileo::Int
#        rdcode::String
#        p0th::T
#        ifprmesh::Bool
#        glnum::Vector{Int}
#        x::Array{T}
#        y::Array{T}
#        z::Array{T}
#        u::Array{T}
#        v::Array{T}
#        w::Array{T}
#        p::Array{T}
#        t::Array{T}
#
#      end
##---------------------------------------------------------------------- 
      mutable struct Re2Hdr 

#        hdr::String
        version::String
        wdsize::Int
        nelgt::Int
        ldim::Int
        nelgv::Int

      end       
#---------------------------------------------------------------------- 

      mutable struct Re2Field{T} #where T <: AbstractFloat

#        wdsize::Int
#        hdr::String
#        version::String
#        nelgt::Int
#        ldim::Int
#        nelgv::Int
        hdr::Re2Hdr
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
     mutable struct Ma2Hdr

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
      mutable struct Ma2Field

#       .ma2 data
#        hdr::String
        hdr::Ma2Hdr        
        pmap::Vector{Int}
        vmap::Array{Int,2}
      end       

#----------------------------------------------------------------------





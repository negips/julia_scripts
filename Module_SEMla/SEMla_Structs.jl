#     Add Definition of new structures here
#---------------------------------------------------------------------- 
"""
      struct SEMla_Config

      Has fields:

      dims::Int             # Physical Dimension of the problem

"""
struct SEMla_Config

  dims::Int

end
#----------------------------------------------------------------------       
"""
      mutable struct TensorField{T <: Number,N} <: AbstractTensorField{T,N}

      Has fields:

      n::Vector{Int}        # Points in the different dimensions
      nel::Int              # No of Tensor elements
      tfield::Array{T}      # Tensor Field

"""
mutable struct TensorField{T <: Number,N} <: AbstractTensorField{T,N}
  n::Vector{Int}
  nel::Int
  tfield::Array{T,N}
end
#---------------------------------------------------------------------- 
"""
      mutable struct IsoTensorField{T <: Number,N} <: AbstractIsoTensorField{T}

      Has fields:

      n::Int                # Constant no of points for all dimensions
      nel::Int              # No of Tensor elements
      tfield::Array{T}      # Tensor Field

"""
mutable struct IsoTensorField{T <: Number, N} <: AbstractIsoTensorField{T,N}
  n::Int
  nel::Int
  tfield::Array{T,N}
end
#---------------------------------------------------------------------- 
@doc raw"""
      mutable struct NTensorFields{T <: Number,N} <: AbstractTensorField{T}

      Has fields:

      n::Vector{Int}              # Points in the different dimensions
      nel::Int                    # No of Tensor elements
      ntensors::Int               # No of Tensor Fields
      TFields::Vector{Array{T}}   # Vector of TensorFields 

"""
mutable struct NTensorFields{T <: Number,N} <: AbstractTensorField{T,N}

  n::Vector{Int}
  nel::Int
  ntensors::Int
  TFields::Vector{TensorField{T,N}}
end
#---------------------------------------------------------------------- 

@doc raw"""
      struct LocalVertexMap
      
      Has fields:

      unique_nvert::{Int}                 - No of Unique vertices
      CI::Vector{CI}                      - Indices of unique vertices sorted in order of their appearance in the array
      reps::Vector{Int}                   - No of repetitions of each vertex
      map2u::AbstractArray{Int}           - Map from original array to unique vector
      Vmap::AbstractArray{Int}            - Vertex Numbering 
                                          - CI: CartesianIndex

"""
struct LocalVertexMap

  unique_nvert::Int
  CI::AbstractVector{CartesianIndex}
  reps::AbstractVector{Int}
  map2u::AbstractArray{Int}
  vmap::AbstractArray{Int}
 
end
#----------------------------------------------------------------------
@doc raw"""
      struct GlobalVertexMap

      Global Vertex Mapping
      
      Has fields:

      npids::Int                          - Total no. of processes with shared edges
      pids::Vector{Int}                   - Process ids
      nedges::Vector{Int}                 - No of edges shared with each process
      map2u::AbstractMatrix{Int}          - Map shared edge (of each process) to local unique edge
      map2p::AbstractMatrix{Int}          - Map local unique edge to shared edges (of each process)
      LocalMap::LocalVertexMap            - Corresponding Local Map

"""
struct GlobalVertexMap
  npids::Int
  pids::AbstractVector{Int}
  nvertex::AbstractVector{Int}
  map2u::AbstractMatrix{Int}
  map2p::AbstractMatrix{Int}
  LocalMap::LocalVertexMap
end
#---------------------------------------------------------------------- 
@doc raw"""
      struct LocalEdgeMap
      
      Has fields:

      unique_nedge::{Int}                       - No of Unique vertices
      unique_edges::Vector{Tuple{Int,Int}}      - Tuple containing (egdeno,elno) of unique edges
     
      reps::Vector{Int}                         - Repetitions of each unique edge
      CI::Vector{CartesianIndex}                - Cartisian Index of the (first) unique elements in the original array
      edgemap::Matrix{Int}                      - Map edges to vector of local unique edges. (Local Edge Number)
      orient::Matrix{Int}                       - Orientation w.r.t. the unique edges

"""
struct LocalEdgeMap
  unique_nedge::Int
  unique_edges::AbstractVector{Tuple{Int,Int}}       
  reps::AbstractVector{Int}
  CI::AbstractVector{CartesianIndex}
  edgemap::AbstractMatrix{Int}
  orient::AbstractMatrix{Int}
 
end
#---------------------------------------------------------------------- 

@doc raw"""
      struct GlobalVertexMap1S

      Global Vertex Mapping for One-Sided communication (1S)
      
      Has fields:

      shared_total_np::Int                    - No of Processes with which data is shared
      shared_nvert::Int                       - No of shared vertices
      vertex_total_np::Vector{Int}            - No of shared vertices (local)
      shared_CI::Vector{CI}                   - Cartesian Indices of shared vertices in original array 
      shared_pid::AbstractMatrix{Int}         - Process ids associated with shared vertices 
      shared_indx::AbstractMatrix{Int}        - Index of shared vertices in a different process id. (one-sided communication) 
      Vmap::AbstractArray{Int}                - Vertex Numbering                                             
                                              - CI: CartesianIndex

"""
struct GlobalVertexMap1S
  shared_total_np::Int
  shared_nvert::Int
  vertex_total_np::Vector{Int}
  shared_CI::Vector{CartesianIndex}
  shared_pid::AbstractMatrix{Int}
  shared_indx::AbstractMatrix{Int}
  Vmap::AbstractArray{Int}
end
#---------------------------------------------------------------------- 

@doc raw"""
      struct GlobalCommMap

      Global Communication Map
      
      Has fields:

      npids::Int                    - No of Processes with which data is shared
      pids::Vector{Int}             - Process Ids
      nverts::Vector{Int}           - No of shared vertices with each process.
      CI::AbstractMatrix{CI}        - Cartesian Indices of shared vertices in original array 
      indx::AbstractMatrix{Int}     - Index of shared vertices in a different process id. (one-sided communication) 
                                    - CI: CartesianIndex


"""
struct GlobalCommMap
  npids::Int
  pids::Vector{Int}
  nverts::Vector{Int}
  CI::AbstractMatrix{CartesianIndex}
  indx::AbstractMatrix{Int}
  Vmap::AbstractArray{Int}
end
#---------------------------------------------------------------------- 
@doc raw"""
      struct GlobalEdgeMap

      Global Edge Mapping
      
      Has fields:

      npids::Int                          - Total no. of processes with shared edges
      pids::Vector{Int}                   - Process ids
      nedges::Vector{Int}                 - No of edges shared with each process
      map2u::AbstractMatrix{Int}          - Map shared edge (of each process) to local unique edge
      map2p::AbstractMatrix{Int}          - Map local unique edge to shared edges (of each process) 
      orient::AbstractMatrix{Int}         - Orientation of shared edge w.r.t unique edge

      We have two distinct maps because order of data in the send and receives is important.

"""
struct GlobalEdgeMap
  npids::Int
  pids::Vector{Int}
  nedges::Vector{Int}
  map2u::AbstractMatrix{Int}
  map2p::AbstractMatrix{Int}
  orient::AbstractMatrix{Int}
end
#---------------------------------------------------------------------- 






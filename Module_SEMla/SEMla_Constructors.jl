#     Add Definition of Constructors here
#---------------------------------------------------------------------- 
"""
      function TensorField(field::Array) where T <: Number
        
      Constructor for the TensorField struct. 

"""
function TensorField(field::Array{T,N}) where {T <: Number,N}

   s    = collect(Int,size(field))
   
   l    = length(s)
   nel  = s[l]
   n    = s[1:l-1]

   return TensorField{T,N}(n,nel,field)
 end        

#---------------------------------------------------------------------- 
"""
      function IsoTensorField(field::Array) where T <: Number
        
      Constructor for the IsoTensorField struct. 

"""
function IsoTensorField(field::Array{T,N}) where {T <: Number, N}

   s    = collect(Int,size(field))
   n    = s[1]
   l    = length(s)
   dims = l-1
   for i in s[1:dims] @assert i == n "Unequal Tensor dimensions $i, $p" end
  
   nel  = s[l]

   return IsoTensorField{T,N}(n,nel,field)
end        

#---------------------------------------------------------------------- 
"""
      function NTensorFields(f::Vararg{Array{T,N}}) where {T<:Number,N}
          
      Constructor for the NTensorFields struct. 

"""
function NTensorFields(f::Vararg{Array{T,N}}) where {T<:Number,N}

  ntensors  = length(f)

  s0        = size(f[1]) 
  for n in 1:ntensors
    si = size(f[n])
    @assert  si == s0 "Unequal Tensor sizes $si ≠ $s0"
  end  
  
  s         = collect(Int,s0)
  l         = length(s)
  nel       = s[l]
  n         = s[1:l-1]

  tfields   = Vector{TensorField{T,N}}(undef,ntensors)
  for m in 1:ntensors
    tfields[m] = TensorField(f[m])
  end  

  return NTensorFields{T,N}(n,nel,ntensors,tfields)
end        

#---------------------------------------------------------------------- 
"""
      function NTensorFields(f::Vector{Array{T,N}}) where {T<:Number,N}
          
      Constructor for the NTensorFields struct. 

"""
function NTensorFields(f::Vector{Array{T,N}}) where {T<:Number,N}

  ntensors  = length(f)

  s0        = size(f[1]) 
  for n in 1:ntensors
    si = size(f[n])
    @assert  si == s0 "Unequal Tensor sizes $si ≠ $s0"
  end  
  
  s         = collect(Int,s0)
  l         = length(s)
  nel       = s[l]
  n         = s[1:l-1]

  tfields   = Vector{TensorField{T,N}}(undef,ntensors)
  for m in 1:ntensors
    tfields[m] = TensorField(f[m])
  end  

  return NTensorFields{T,N}(n,nel,ntensors,tfields)
end        

#---------------------------------------------------------------------- 
"""
      function NTensorFields(f::Vector{TensorField{T,N}}) where{T<:Number,N}
          
      Constructor for the NTensorFields struct. 

"""
function NTensorFields(f::Vector{TensorField{T,N}}) where {T<:Number,N}

  ntensors  = length(f)

  nel0      = f[1].nel
  n0        = f[1].n
  for n in 1:ntensors
    neli = f[n].nel
    ni   = f[n].n
    @assert  neli == nel0 "Unequal TensorField nel $neli ≠ $nel0"
    @assert  ni   == n0   "Unequal TensorField n $ni ≠ $n0"
  end  

  tfields = f

  return NTensorFields{T,N}(n0,nel0,ntensors,tfields)
end        

#---------------------------------------------------------------------- 
"""

      LocalVertexMap(vmap::Array{Int})

      Constructor for LocalVertexMap struct.
      Get the local unique vertices (and their mapping) from the vertex map (vmap).

"""
function LocalVertexMap(vmap::Array)


  nunq,indu_s,reps_s,rind_s,rindu_s = localmapping(vmap)

  l                           = length(vmap)
  unique_CI                   = Vector{CartesianIndex}(undef,nunq)
  repeated_CI                 = Vector{CartesianIndex}(undef,l)
  unique_in_repeated_I        = Vector{Int}(undef,nunq)

  CI = CartesianIndices(vmap)
  for i in 1:nunq
    unique_CI[i]              = CI[indu_s[i]]
    unique_in_repeated_I[i]   = rindu_s[i]
  end  

  for i in 1:l
    repeated_CI[i]  = CI[rind_s[i]]
  end  

  localmap = LocalVertexMap(nunq,unique_CI,reps_s,repeated_CI,unique_in_repeated_I)    

#  localmap = LocalVertexMap(nunq,indu_s,reps_s,rind_s,rindu_s)    

  return localmap
end
##---------------------------------------------------------------------- 
"""

      LocalVertexMap()

      Constructor for empty LocalVertexMap object.

"""
function LocalVertexMap()

  unique_nvert    = 0
  CI              = Vector{CartesianIndex}(undef,0)
  reps            = Vector{Int}(undef,0)
  map2u           = Matrix{Int}(undef,0,0)
  vmap            = Matrix{Int}(undef,0,0)

  localmap        = LocalVertexMap(unique_nvert,CI,reps,map2u,vmap)    

  return localmap
end
##---------------------------------------------------------------------- 

"""

      GlobalVertexMap1S()

      Constructor for empty GlobalVertexMap1s.
      Needed when we have no shared vertices.
      for np = 1 for example.

"""
function GlobalVertexMap1S()

  total_np = 0
  nvert    = 0
  v_np     = Vector{Int}(undef,0)
  CI       = Vector{CartesianIndex}(undef,0)
  pid      = Matrix{Int}(undef,0,0)
  indx     = Matrix{Int}(undef,0,0)

  return GlobalVertexMap1S(total_np,nvert,v_np,CI,pid,indx)
end
##---------------------------------------------------------------------- 
"""

      GlobalVertexMap()

      Constructor for empty GlobalVertexMap1s.
      Needed when we have no shared vertices.
      for np = 1 for example.

"""
function GlobalVertexMap()

  npids    = 0
  pids     = Vector{Int}(undef,0)
  nvertex  = Vector{Int}(undef,0)
  map2u    = Matrix{Int}(undef,0,0)
  map2p    = Matrix{Int}(undef,0,0)
  LocalMap = LocalVertexMap()

  return GlobalVertexMap(npids,pids,nvertex,map2u,map2p,LocalMap)
end
##---------------------------------------------------------------------- 

"""

      GlobalCommMap()

      Constructor for empty GlobalCommMap.
      Needed when we have no shared vertices.
      for np = 1 for example.

"""
function GlobalCommMap()

  npids    = 0
  pids     = Vector{Int}(undef,0)
  nverts   = Vector{Int}(undef,0)
  CI       = Matrix{CartesianIndex}(undef,0,0)
  indx     = Matrix{Int}(undef,0,0)

  return GlobalCommMap(npids,pids,nverts,CI,indx)
end
##---------------------------------------------------------------------- 















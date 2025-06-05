# Storing old unused functions here

"""

      buildlocalmap1(vmap::Array{Int})

      Build a LocalVertexMap object from vmap
      Get the local unique vertices (and their mapping) from the vertex map.

"""
function buildlocalmap1(vmap::Array)


  nunq,indu_s,reps_s,rind_s,rindu_s = localmapping1(vmap)

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

  lmap = LocalVertexMap(nunq,unique_CI,reps_s,repeated_CI,unique_in_repeated_I)    

  return lmap
end
##---------------------------------------------------------------------- 
"""

      localmapping1(vmap::Array{Int})

      Generate the local unique vertices (and their mapping) from the vertex map (vmap).

"""
function localmapping1(vmap::AbstractArray{Int})

  s0        = ndims(vmap)
  dims      = s0 - 1
  l         = length(vmap)

  v0        = vmap[:]        
  ind       = sortperm(v0)
  indi      = invperm(ind)    # Inverted permutation

  v1        = Base.copy(v0[ind])   # And v0 == v1[indi]
  T         = eltype(v1)
  reps      = zeros(T,l)      # No of repetitions of vertex maps
  indu      = zeros(T,l)      # Indices of unique points
  rind      = zeros(T,l)      # Indices of repeated vertices
  rindu     = zeros(T,l)      # Index of indu elements in rind array

  # Count repetitions and total no of unique vertices
  n         = 1
  num       = v1[1]
  st        = 1
  en        = 1
  nunq      = 1               # No of (local) unique vertices
  for i in 2:l
    num2 = v1[i]
    if num2 == num
      # Repetition
      n  = n+1
      en = i
    else
      # New vertex found
      si   = st:en            # Sorted Indices
      oi   = indi[si]         
      ni   = ind[si]

      ri   = view(rind, si)
      Base.copy!(ri,ni)
      j    = minimum(ni)      # take the first index in the original array
      indu[nunq] = j
      reps[nunq] = n
      rindu[nunq]= st         # Pointer of indu element in rind array
     
      st   = i
      en   = i
      n    = 1
      num  = num2
      nunq = nunq+1
    end

    if i == l
      si   = st:en
      oi   = indi[si]
      ni   = ind[si]
      ri   = view(rind, si)
      Base.copy!(ri,ni)
      j    = minimum(ni)      # take the first index in the original array
      indu[nunq] = j
      reps[nunq] = n
      rindu[nunq]= st         # Pointer of indu element in rind array
    end
  end  

  # In principle this sorting is not necessary
  # But it should reduce jumps in memory access
  i4        = sortperm(indu[1:nunq])
  i5        = invperm(i4) 

  indu_s    = Base.copy(indu[i4])
  reps_s    = Base.copy(reps[i4])
  rindu_s   = Base.copy(rindu[i4])

  reps_c    = cumsum(reps_s)
  for i in 1:nunq
    reps_c[i] = reps_c[i] - reps_s[i]
  end  

  reps_cu   = cumsum(reps)
  for i in 1:nunq
    reps_cu[i] = reps_cu[i] - reps[i]
  end  

  rind_s    = Base.copy(rind)
  for i in 1:nunq
    j   = i5[i]         # j - New position of the old ith element
    rj  = reps_c[j]     # Offset due to previous elements
    nr  = reps_s[j]     # No of repeated indices of jth unique element
    st  = rj+1          # start position in the new sorted vector

    ri  = reps_cu[i]    # Offset due to previous elements (In original vector)
    sti = ri + 1

    Base.copyto!(rind_s,st,rind,sti,nr)
    
    rindu_s[j] = st

  end  

  # create the structure
  # I do structure creation within the constructor.

  return nunq,indu_s,reps_s,rind_s,rindu_s
end

#----------------------------------------------------------------------

"""

      buildglobalmap1(vmap::Array{Int},localmap::LocalVertexMap,comm::MPI.Comm)

      Generate a global mapping from the Vertex map (vmap) and the local mapping (localmap).

      Input:
      vmap        - Vertex map
      localmap    - LocalVertexMap object
      comm        - MPI Communicator

      Output:
      globalmap   - GlobalVertexMap1s object

"""
function buildglobalmap1(vmap::AbstractArray{Int},localmap::LocalVertexMap,comm::MPI.Comm)

  local np    = MPI.Comm_size(comm)
  local prank = MPI.Comm_rank(comm)

  l         = ndims(vmap)
  dims      = l - 1
  repsmax   = 2^dims    # Max no of repetitions possible.

  # Declare initial array for global vertices
  gvertex   = fill(-1,localmap.unique_nvert)
  CI        = Vector{CartesianIndex}(undef,localmap.unique_nvert)

  nloc1 = 0
  for i in eachindex(localmap.unique_reps)
    if localmap.unique_reps[i]<repsmax
      nloc1         += 1
      ci             = localmap.unique_CI[i]
      CI[nloc1]      = ci
      gvertex[nloc1] = vmap[ci]
    end
  end  

  # Total number of globally shared vertices  
  nglob     = MPI.Allreduce(nloc1, +, comm)            # Global no of points
  nlocmax   = MPI.Allreduce(nloc1, MPI.MAX, comm)      # maximum local points

  # Return an empty globalmap if there are no vertices to share.
  # Full periodic domain on single processor
  if nlocmax == 0
    return GlobalVertexMap1S()
  end  

  # Assume we know this for now. Declared Constant Possibly?
  N1    = 10
  # Usual size of a local array
  npts  = length(vmap)*(N1^dims)     

  # Assume we have a scratch space that is 3 times the size
  npts2 = npts*3 

  nvmax       = 2^dims - 1

  if nlocmax*np < npts2
    # All points fit within one array.
    if prank == 0
      println("Everything fits in one array, $(nlocmax*np), $npts2")
    end

    wka = fill(-1,npts2)

    s_np,ns,s_pid,s_indx = localtoglobalmapping(gvertex,nloc1,dims,wka,comm)

    # Thow out vertices which have no global connections
    ns2 = 0
    for k in 1:nloc1
      if ns[k] > 0
        ns2 += 1
      end
    end

    ns2max  = MPI.Allreduce(ns2, MPI.MAX, comm)      # maximum local points

    # If we have no shared vertices, return an empty global map
    if ns2max == 0
      gmap = GlobalVertexMap1S()
    else  

      CI2     = Vector{CartesianIndex}(undef,ns2)
      # gvert2  = view(gvertex,1:ns2max)            # Reuse gvertex
      gvert2  = fill(-1,ns2max)

      i = 0
      for k in 1:nloc1
        if ns[k] > 0
          i += 1
          gvert2[i] = gvertex[k]
          CI2[i]    = CI[k]
        end
      end

      for k in nloc1+1:localmap.unique_nvert
        gvertex[k] = -1
      end  

      # Generate the map again from the trimmed data set
      nloc2 = ns2
   
      nps_tot,shared_np,shared_pid,shared_indx = localtoglobalmapping(gvert2,nloc2,dims,wka,comm)

      # create object
      gmap = GlobalVertexMap1S(nps_tot,nloc2,shared_np,CI2,shared_pid,shared_indx,vmap)

      # comm_map  = globalcomm_map(globalmap)

    end           # ns2max == 0 
 
  else
    # (nlocmax*np) > npts2
    #
    # To be implemented. Probably need to repeat the above process in chunks

  end       # (nlocmax*np) < npts2


  MPI.Barrier(comm)

  println("$(gmap.shared_total_np) processes shared by Rank $prank\n")

  return gmap
end
#---------------------------------------------------------------------- 
"""

      localtoglobalmapping(gvert::Vector{Int},nloc::Int,dims::Int,wk::Vector{Int},comm::MPI.Comm)

      Generate a global mapping from the vertex numbers gvert.

      Input:
      nloc  - No. of local vertices
      dims  - dimenson of physical space
      wk    - work array
      comm  - MPI Communicator

      Output:
      nproc_share::Int              - No. of Process sharing vertices with this prank    
      nshared::Vector{Int}          - No. of Processes sharing each vertex
      shared_pid::Matrix{Int}       - Process ranks (multiple) of earch shared vertex
      shared_indx::Matrix{Int}      - Index in shared array of the corresponding process

"""
function localtoglobalmapping(gvert::AbstractVector{Int},nloc::Int,dims::Int,wk::Vector{Int},comm::MPI.Comm)

  # All points fit within one array.
  
  local prank = MPI.Comm_rank(comm)
  local np    = MPI.Comm_size(comm)

  # maximum local points
  nlocmax   = MPI.Allreduce(nloc, MPI.MAX, comm)

  t1        = nlocmax*np
  t2        = length(wk)
  @assert t2>t1 "Insufficient work array length: $t1, $t2"

  fill!(wk,-1)
  # offs     = prank*nlocmax+1

  # println("$prank: Array: $(nlocmax*np),  el: $(nlocmax)")

  vw       = view(gvert,1:nlocmax)
  vw_all   = view(wk,1:nlocmax*np)
  # gvert -> wk
  MPI.Allgather!(vw,vw_all,comm)    

  # Allocate
  nvmax       = 2^dims - 1
  # Process id of each shared point 
  shared_pid  = fill(-1,nloc,nvmax)

  # Index in the shared vector
  shared_indx = fill(-1,nloc,nvmax)       

  # No of Processes sharing each vertex
  nshared     = fill(0,nloc)


  # No of Processors with which prank has shared vertices.
  nproc_share = 0

  for pid in 0:np-1
   
    if pid == prank
      continue  # skip this loop
    end
    offs  = pid*nlocmax
    # no of shared vertices with each pid
    vs    = 0   
    for i in 1:nlocmax
      j     = offs+i
      vi    = wk[j]

      if vi<0
        break
      end  

      for k in 1:nloc
        if vi == gvert[k]
          ns                = nshared[k]+1
          nshared[k]        = ns
          shared_pid[k,ns]  = pid
          shared_indx[k,ns] = i
          vs = vs + 1
        end
      end       # k in 1:nloc
    end         # i in 1:nlocmax
    
    # Found shared points?
    if vs > 0
      nproc_share = nproc_share+1
    end  
  end     # pid in 0:np-1 
 
  return nproc_share,nshared,shared_pid,shared_indx
end
#----------------------------------------------------------------------

@doc raw"""

      globalcomm_map(glmap::GlobalVertexMap1S)

      Build Communication map based on the Global Vertex Map

"""
function globalcomm_map(globalmap::GlobalVertexMap1S)

  # Rearrange the global vertex map by process id for easier communication
  # There's no new information. Just easier access.

  nps_tot = globalmap.shared_total_np

  if nps_tot == 0
    return GlobalCommMap()
  end  


  # Unique pids
  pidu   = fill(-1,nps_tot)

  # Total vertices shared with each pid
  pverts = fill(0,nps_tot) 

  npid        = 0
  for i in eachindex(globalmap.shared_pid)
    p   = globalmap.shared_pid[i]
    if p>-1
      j   = findfirst(x -> x==p, pidu)
      if isnothing(j)
        npid    += 1
        j        = npid
        pidu[j]  = p
      end
      pverts[j] += 1
    end  
  end  
  @assert npid == nps_tot "Something went wrong in getting unique pids"

  pvmax     = maximum(pverts)
  local CI  = Matrix{CartesianIndex}(undef,pvmax,nps_tot) 
  indx      = Matrix{Int}(undef,pvmax,nps_tot) 

  lverts    = fill(0,nps_tot)
  for ci in CartesianIndices(globalmap.shared_pid)
    p       = globalmap.shared_pid[ci]
    tup     = Tuple(ci)
    k       = tup[1]
    if p < 0
      continue
    end
    j          = findfirst(x -> x==p, pidu)
    lverts[j] += 1
    l          = lverts[j]
    CI[l,j]    = globalmap.shared_CI[k]
    indx[l,j]  = globalmap.shared_indx[ci]
  end  

  @assert pverts == lverts "Something went wrong in counting vertices"

  comm_map = GlobalCommMap(npid,pidu,pverts,CI,indx,globalmap.Vmap)

  return comm_map 
end

#---------------------------------------------------------------------- 

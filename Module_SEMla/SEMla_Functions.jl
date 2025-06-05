#     Non-Constructor Function Definitions
#---------------------------------------------------------------------- 
#----------------------------------------------------------------------       
"""

      buildlocalvertexmap(vmap::Array{Int})

      Build a LocalVertexMap object from vmap
      Get the local unique vertices (and their mapping) from the vertex map.

"""
function buildlocalvertexmap(vmap::Array)


  vmapu,reps,CI,vertmap = localuniquemap(vmap)
  nunq = length(vmapu)

  lmap = LocalVertexMap(nunq,CI,reps,vertmap,vmap)    

  return lmap
end
##---------------------------------------------------------------------- 
"""

      buildglobalvertexmap(localmap::LocalVertexMap,comm::MPI.Comm)

      Generate a global mapping from the LocalVertexMap object (localmap).

      Input:
      localmap    - LocalVertexMap object (already contains vmap)
      comm        - MPI Communicator

      Output:
      globalmap   - GlobalVertexMap object

"""
function buildglobalvertexmap(localmap::LocalVertexMap,comm::MPI.Comm)

  local np    = MPI.Comm_size(comm)
  local prank = MPI.Comm_rank(comm)

  vmap        = localmap.vmap
  l           = ndims(vmap)
  dims        = l - 1
  repsmax     = 2^dims    # Max no of repetitions possible.
  
  # In unstructured meshes we can have more.
  repsmax     = 100


  numax   = MPI.Allreduce(localmap.unique_nvert, MPI.MAX, comm)      # maximum local points

  # Declare initial array for global vertices
  gvertex     = fill(-1,numax)
  sharedtounq = fill(-1,numax)

  nloc1 = 0
  for i in 1:localmap.unique_nvert
    if localmap.reps[i]<repsmax
      nloc1             += 1
      ci                 = localmap.CI[i]
      gvertex[nloc1]     = vmap[ci]
      sharedtounq[nloc1] = i
    end
  end  

  # Total number of globally shared vertices  
  nglob     = MPI.Allreduce(nloc1, +, comm)            # Global no of points
  nlocmax   = MPI.Allreduce(nloc1, MPI.MAX, comm)      # maximum local points

  # Return an empty globalmap if there are no vertices to share.
  # Full periodic domain on single processor
  if nlocmax == 0
    return GlobalVertexMap()
  end  

  # Assume we know this for now. Declared Constant Possibly?
  # Or use a sub module for scratch arrays?
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

    # work array
    wka = fill(-1,npts2)

    npids,pids,eshare,maptou,maptop = localuniquetoglobal(gvertex,nloc1,wka,comm)
    for ci in CartesianIndices(maptou)
      ue1 = maptou[ci]
      if ue1>-1
        ue2 = sharedtounq[ue1]
        maptou[ci] = ue2
      end  

      ue1 = maptop[ci]
      if ue1>-1
        ue2 = sharedtounq[ue1]
        maptop[ci] = ue2
      end  
    end  

    gmap = GlobalVertexMap(npids,pids,eshare,maptou,maptop,localmap)
 
  else
    # (nlocmax*np) > npts2
    #
    # To be implemented. Probably need to repeat the above process in chunks

    error("Multiprocessor global vertex mapping not yet implemented.")

  end       # (nlocmax*np) < npts2


  MPI.Barrier(comm)

  println("$(gmap.npids) processes shared by Rank $prank\n")

  return gmap
end
#---------------------------------------------------------------------- 
@doc raw"""

      buildlocaledgemap(vmap::AbstractArray{Int})

      Get local unique edges (by vertex numbers),
      repetitions and mapping.

"""
function buildlocaledgemap(vmap::AbstractArray{Int})

  local sedges,edgeu,edgemap,reps,CI,ledgemap,edgenos

  nd     = ndims(vmap)
  dims   = nd - 1

  sedges = sorted_edgevertices(vmap)

  # unique edges, repetitions,mapping
  edgeu,reps,CI,edgemap = localuniquemap(sedges)
  neu       = length(edgeu)

  edgenos   = Vector{Tuple{Int,Int}}(undef,neu)
  for i in 1:neu
    ci         = CI[i]
    tup        = Tuple(ci)
    edgenos[i] = tup
  end

  # Edge Orientations
  edgeorient = fill(0,size(edgemap))
  for ci in CartesianIndices(edgeorient)
    i       = edgemap[ci]     # index in unique array
    tupu    = edgenos[i]
    edu     = tupu[1]
    el      = tupu[2]
    e1,e2   = edgeindices(edu,dims)
    c1      = CartesianIndex(e1)
    c2      = CartesianIndex(e2)
    # Vertices of the unique edge
    vertu   = [vmap[c1,el]; vmap[c2,el]]

    tup     = Tuple(ci)
    edno    = tup[1]
    el2     = tup[2]
    e3,e4   = edgeindices(edno,dims)
    c3      = CartesianIndex(e3)
    c4      = CartesianIndex(e4)
    vert    = [vmap[c3,el2]; vmap[c4,el2]]

    if vertu[1] == vert[1] && vertu[2] == vert[2]
      edgeorient[ci] = 1
    elseif vertu[1] == vert[2] && vertu[2] == vert[1]
      edgeorient[ci] = -1
    else
      error("Something went wrong.")
    end
  end  

  ledgemap  = LocalEdgeMap(neu,edgenos,reps,CI,edgemap,edgeorient)

  return ledgemap
end
#---------------------------------------------------------------------- 
@doc raw"""

      sorted_edgevertices(vmap::Array{Int})

      Get array of edge vertex tuples. 
      Each tuple is sorted by vertex number.

"""
function sorted_edgevertices(vmap::Array{Int})

  s         = size(vmap)
  l         = ndims(vmap)
  nel       = s[l]
  dims      = l - 1

  if dims == 1
    nedges = 0
  elseif dims == 2
    nedges = 4
  elseif dims == 3
    nedges = 12
  end

  s_edge_verts = Matrix{Tuple{Int,Int}}(undef,nedges,nel)

  for e in 1:nel
    for i in 1:nedges
      tup   = edgeindices(i,dims)
      ci1   = CartesianIndex(tup[1])
      ci2   = CartesianIndex(tup[2])
      ve    = [vmap[ci1,e]; vmap[ci2,e]]
      s_ve  = sort(ve)

      s_edge_verts[i,e] = tuple(s_ve[1], s_ve[2])
    end
  end
   
  return s_edge_verts
end
##----------------------------------------------------------------------
"""

      buildglobaledgemap(ledgemap::LocalEdgeMap,gvmap::GlobalVertexMap,dims::Int,comm::MPI.Comm)

      Generate a global mapping from the vertex numbers gvert.

      Input:
      nloc  - No. of elements in vector
      dims  - dimenson of physical space
      wk    - work array
      comm  - MPI Communicator

      Output:
      nproc_share::Int              - No. of Process sharing vertices with this prank    
      nshared::Vector{Int}          - No. of Processes sharing each vertex
      shared_pid::Matrix{Int}       - Process ranks (multiple) of earch shared vertex
      shared_indx::Matrix{Int}      - Index in shared array of the corresponding process

"""
function buildglobaledgemap(ledgemap::LocalEdgeMap,gvmap::GlobalVertexMap,dims::Int,comm::MPI.Comm)
  
  local prank = MPI.Comm_rank(comm)
  local np    = MPI.Comm_size(comm)
  local nloc, nlocmax 

  if dims == 1
    emax = 0
  elseif dims == 2
    emax = 2
  elseif dims == 3
    emax = 4
  end

  # In unstructured meshes we can have more
  emax = 100

  # Get index of elements which could possibly be shared.
  # Less work in the next loop
  ne = 0
  for i in 1:ledgemap.unique_nedge
    ereps = ledgemap.reps[i]
    if ereps < emax
      ne += 1
    end
  end
  neshared    = ne

  nemax       = MPI.Allreduce(ne, MPI.MAX, comm)

  # Set up memory to receive the edges
  snp         = gvmap.npids
  sharededges = Matrix{Tuple{Int,Int}}(undef,nemax,snp) 
  rstat       = Vector{Any}(undef,snp)

  for i in 1:snp
    p         = gvmap.pids[i]
    vw        = view(sharededges,:,i)
    rstat[i]  = MPI.Irecv!(vw,comm; source=p)
  end  

  localedges   = Vector{Tuple{Int,Int}}(undef,nemax)

  vmap         = gvmap.LocalMap.vmap
  ne           = 0
  for j in 1:ledgemap.unique_nedge
    ereps = ledgemap.reps[j]
    if ereps < emax
      ne                += 1
      etup               = ledgemap.unique_edges[j]
      edgno,el           = collect(Int,etup)
      e1,e2              = edgeindices(edgno,dims)
      c1                 = CartesianIndex(e1)
      c2                 = CartesianIndex(e2)
      localedges[ne]     = tuple(vmap[c1,el], vmap[c2,el])
    end
  end  

  # send buffer is just a view
  vw  = view(localedges, 1:ne)

  # Send Vertex Tuples
  for i in 1:snp
    p  = gvmap.pids[i]
    tag = i
    MPI.Isend(vw,comm;dest=p,tag=i)
  end

  # Wait for the Receives
  recv_count = fill(0,snp)
  for k in 1:snp
    stat = MPI.Wait!(rstat[k])
    recv_count[k] = MPI.Get_count(stat,Tuple{Int,Int})
    # println("Received ",recv_count[k]," tuples from ",gvmap.pids[k], " at rank ", prank)
  end  

  maptou          = fill(-1,nemax,snp)
  maptop          = fill(-1,nemax,snp)
  edgeorient      = fill(0,size(maptou))
  sharedpids      = fill(-1,snp)
  shared_nedges   = fill(0,snp)

  vec_edgemap     = fill(-1,nemax)
  vec_orient      = fill(-1,nemax)

  sp = 0    # no. of shared processors
  for i in 1:snp
    l = 0
    for j in 1:recv_count[i]
      tupg = sharededges[j,i]
      v1   = tupg[1]
      v2   = tupg[2]
      for k in 1:ledgemap.unique_nedge
        # Compare with local edge vertices
        etup      = ledgemap.unique_edges[k]
        edgno,el  = collect(Int,etup)
        e1,e2     = edgeindices(edgno,dims)
        c1        = CartesianIndex(e1)
        c2        = CartesianIndex(e2)
        # Local unique edge vertices: w1, w2
        w1        = vmap[c1,el]
        w2        = vmap[c2,el]
        if v1==w1 && v2==w2
          l                   +=  1
          vec_orient[l]        =  1
          vec_edgemap[l]       =  k
        elseif  v1==w2 && v2==w1
          l                   +=  1
          vec_orient[l]        = -1
          vec_edgemap[l]       =  k
        end
      end
    end     # j in 1:recv_count[i]

    #if l>0
      # Common edges were found
      # Increase processor count for shared edges
      sp    += 1
      sharedpids[i]     = gvmap.pids[i]
      shared_nedges[i]  = l

      vw = view(maptou,:,i)
      copyto!(vw,1,vec_edgemap,1,l)

      vw2 = view(edgeorient,:,i)
      copyto!(vw2,1,vec_orient,1,l)
    #end     # l > 0
  end       # i = 1:gvmap.npids


  # Build the reverse map.
  # Unique edges to processes
  for i in 1:snp
    nf = 0
    for j in 1:ledgemap.unique_nedge
      found = false
      for k in 1:shared_nedges[i]
        if j == maptou[k,i]
          found = true
          break
        end
      end
      if found
        nf += 1
        maptop[nf,i] = j
      end
    end
  end

  gemap = GlobalEdgeMap(sp,sharedpids,shared_nedges,maptou,maptop,edgeorient)     

  return gemap
end
#----------------------------------------------------------------------

@doc """
      function localuniquemap(array::AbstractArray{T}) where {T}

      Get the unique elements of an array.

      Outputs:
      
      array_u::Vector{T}            - Unique elements of array.
      reps::Vector{Int}             - No of repetitions of each unique element.
      CI::Vector{CartesianIndex}    - Cartesian Index of first occurence of each unique element.
      map::Array{Int}               - Map from array -> array_u

"""
function localuniquemap(array::AbstractArray{T}) where {T}

  local array_u,reps,map,CI

  stup   = size(array)
  l      = length(array)
  
  # Unique elements of the array
  array_u   = similar(array,l) # Vector{T}(undef,l)

  reps      = similar(array,Int,l)              # Vector{Int}(undef,l)
  map       = similar(array,Int,stup)           # Array{Int}(undef,s)
  CI        = similar(array,CartesianIndex,l)   # Vector{CartesianIndex}(undef,l) 

  # No of unique elements
  neu = 0
  for i in CartesianIndices(array)
    e = array[i]
    newelement = true
    for k in 1:neu
      eu = array_u[k]
      if e == eu
        newelement = false
        map[i]     = k
        reps[k]   += 1
        break
      end
    end

    if newelement
      neu         += 1
      map[i]       = neu
      array_u[neu] = e
      reps[neu]    = 1
      CI[neu]      = i
    end
     
  end  # i in CartesianIndices(array)

  array_u2 = similar(array_u,neu) # Vector{T}(undef,neu)
  Base.copyto!(array_u2,1,array_u,1,neu)

  reps2  = similar(reps,Int,neu)  # Vector{Int}(undef,neu)
  Base.copyto!(reps2,1,reps,1,neu)

  CI2    = similar(CI,neu) # Vector{CartesianIndex}(undef,neu)
  Base.copyto!(CI2,1,CI,1,neu)

  return array_u2,reps2,CI2,map
end
#---------------------------------------------------------------------- 
"""

      localuniquetoglobal(lvector::AbstractVector{T},nloc::Int,wk::Vector{T},comm::MPI.Comm) where {T}

      Generate a global mapping from the local vector: lvector

      Input:
      lvector     - local vector
      nloc        - No. of local vertices
      wk          - work array
      comm        - MPI Communicator

      Output:
      nproc::Int                - No. of Process sharing elements with this prank    
      pids::Vector{Int}         - Processes ids of neighbors
      eshared::Vector{Int}      - No of elements shared with with process 
      maptou::Matrix{Int}       - Map from shared processs elements to local unique element 
      maptop::Matrix{Int}       - Map from local unique elements to shared processes

      Note: Ensure wk array is empty or already initialized.
"""
function localuniquetoglobal(lvector::AbstractVector{T},nloc::Int,wk::AbstractVector{T},comm::MPI.Comm) where {T}

  # All points fit within one array.
  
  local prank = MPI.Comm_rank(comm)
  local np    = MPI.Comm_size(comm)

  # maximum local points
  nlocmax   = MPI.Allreduce(nloc, MPI.MAX, comm)

  t1        = nlocmax*np
  t2        = length(wk)
  @assert t2>t1 "Insufficient work array length: $t1, $t2"

  # For a more general code. Make sure wk is empty
  # fill!(wk,-1)

  vw       = view(lvector,1:nlocmax)
  vw_all   = view(wk,1:nlocmax*np)
  # gvert -> wk
  MPI.Allgather!(vw,vw_all,comm)    

  # Processes id of neighbors
  pshared  = fill(0,nloc)
  # No of elements shared with each process
  nelshared = fill(0,nloc)

  # Total no of neighboring processors
  nproc = 0

  vw2      = view(lvector,1:nloc)

  for pid in 0:np-1
    vs    = 0   
    if pid == prank
      continue  # skip this loop
    end
    offs  = pid*nlocmax
    # no of shared elements with each pid
    for i in 1:nlocmax
      j     = offs+i
      vi    = wk[j]

      k = findfirst(x -> x==vi,vw2)
      if !isnothing(k)
        vs += 1
      end
    end

    if vs > 0
      # Then we do have shared elements
      nproc += 1
      pshared[nproc]   = pid
      nelshared[nproc] = vs
    end  
  end     # pid in 0:np-1 

  emax    = maximum(nelshared[1:nproc])
  maptou  = fill(-1,emax,nproc)
  maptop  = fill(-1,emax,nproc)

  for l in 1:nproc
    pid = pshared[l]
    vs    = 0   

    offs  = pid*nlocmax
    # no of shared elements with each pid
    for i in 1:nlocmax
      j     = offs+i
      vi    = wk[j]

      k     = findfirst(x -> x==vi,vw2)
      if !isnothing(k)
        vs        += 1
        maptou[vs,l] = k
      end
    end

    vs = 0
    vw3     = view(wk,offs+1:offs+nlocmax)
    for i in 1:nloc
      vi    = lvector[i]
      k     = findfirst(x -> x==vi,vw3)
      if !isnothing(k)
        vs        += 1
        maptop[vs,l] = i
      end
    end  
  end

  pids      = fill(-1,nproc)
  eshared   = fill(0,nproc)
  copy!(pids,pshared[1:nproc])
  copy!(eshared,nelshared[1:nproc])

  return nproc,pids,eshared,maptou,maptop
end
#----------------------------------------------------------------------

@doc raw"""

      eltoprank(elno::Int,bsize::Int)

      Get process id for element with number `elno` based on size of elements per processors `bsize`

"""
function eltoprank(elno::Int,bsize::Int)

  pid = floor(Int,elno/bsize)

  return pid 
end

#---------------------------------------------------------------------- 









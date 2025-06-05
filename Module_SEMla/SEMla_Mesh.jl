#     Author:     Prabal Negi
#     Description: Mesh Related functions
"""
      gen_rema2(case::String)

      Read the case.re2 and case.ma2 files  and generate a case.rema2
      file with elements sorted according to the numbering in case.ma2.
      Boundary and Curve data also made consistent.

      Todo: Add option to output sorted.re2 and sorted.ma2

"""
      function gen_rema2(case::String, ifsorted::Bool, comm::MPI.Comm, nid0::Int)

        local map, re2

        prank = MPI.Comm_rank(comm)

        re2file = case*".re2"        # String concatenation
        ma2file = case*".ma2"        # String concatenation

        map         = read_ma2(ma2file,comm,nid0)
        re2         = read_re2(re2file,comm,nid0)

        MPI.Barrier(comm)

        newre2      = Re2Field()
        newmap      = Ma2Field()

        if prank == nid0

          nelre2         = re2.hdr.nelgt
          nelma2         = map.hdr.nel
          @assert nelre2 == nelma2 "Total Elements don't Match. $(nelre2), $(nelma2)"

          if (ifsorted)          
            # Get index key for sorted map          
            key         = sortperm(map.pmap)
            keyinv      = Base.copy(key)          # Inverse key
            for i in 1:nelre2
              j         = key[i]
              keyinv[j] = i
            end
          else
            key         = Int.(LinearIndices(1:nelre2))
            keyinv      = Base.copy(key)
          end  

          newre2      = copy(re2)

          # Convert to Symmetric numbering          
          ldim    = re2.hdr.ldim
          nc      = 2^ldim
          p2s     = zeros(Int64,nc)
          for i in 1:nc
            p2s[i] = preproctosymm(i)
          end  

          # Renumber Coordinates and Boundary conditions
          for i in 1:nelre2
            j              = key[i]
            newre2.xc[:,i] = re2.xc[p2s,j]
            newre2.yc[:,i] = re2.yc[p2s,j]
            if re2.hdr.ldim == 3
              newre2.zc[:,i] = re2.zc[p2s,j]
            else
              newre2.zc[1,1] = re2.zc[1,1]
            end
            newre2.cbl[:,i]  = re2.cbl[:,j]     # Need to change face numbers
            newre2.bl[:,:,i] = re2.bl[:,:,j]    # Might need to change parameters
          end

          # Renumber Curve ieg        
          for i in 1:re2.ncurve
            el                  = re2.curveieg[i]
            newel               = keyinv[el]
            newre2.curveieg[i]  = newel
          end  

          newmap                = copy(map)
          for i in 1:nelre2
            j                   = key[i] 
            newmap.pmap[i]      = map.pmap[j]
            newmap.vmap[:,i]    = map.vmap[:,j]
          end 

          # Generate "case.rema2.h5" file
          gen_rema2_h5(case,newre2,newmap,ifsorted,comm,nid0)

        end       # prank == nid0  

#       For now just returning the new data        
        return newre2,newmap 
      end
#----------------------------------------------------------------------
      function gen_rema2_h5(case::String,re2::Re2Field,map::Ma2Field,ifsorted::Bool,comm::MPI.Comm,nid0::Int)

        prank = MPI.Comm_rank(comm)

        if prank == nid0 
          nelre2       = re2.hdr.nelgt
          nelma2       = map.hdr.nel
          @assert map.hdr.nel == re2.hdr.nelgt "Total Elements don't Match. $(nelma2), $(nelre2)"

#         Write out an hdf5 file
          fname = case*".rema2.h5"
          fid   = h5open(fname, "w")
          # .re2 data        
          g  = create_group(fid,"Re2")
          g1 = create_group(g,"Params") 
          g2 = create_group(g,"Data")        

          # Parameters
          write_dataset(g1,"ndim",re2.hdr.ldim)
          write_dataset(g1,"nelgv",re2.hdr.nelgv)
          write_dataset(g1,"nelgt",re2.hdr.nelgt)
          write_dataset(g1,"wdsize",re2.hdr.wdsize)
          write_dataset(g1,"ncurve",re2.ncurve)

          # Data
          write_dataset(g2,"xc",re2.xc)
          write_dataset(g2,"yc",re2.yc)
          write_dataset(g2,"zc",re2.zc)
          write_dataset(g2,"bcs",re2.cbl)
          write_dataset(g2,"bcparams",re2.bl)
          if (re2.ncurve>0)
            write_dataset(g2,"curveieg",re2.curveieg)
            write_dataset(g2,"curveiside",re2.curveiside)
            write_dataset(g2,"curveparam",re2.curveparam)
            write_dataset(g2,"curvetype",re2.curvetype)
          end  

          # .ma2 
          h  = create_group(fid,"Ma2")
          h1 = create_group(h,"Params")
          h2 = create_group(h,"Data")

          # Parameters
          write_dataset(h1,"d2",map.hdr.d2)
          write_dataset(h1,"depth",map.hdr.depth)
          write_dataset(h1,"nactive",map.hdr.nactive)
          write_dataset(h1,"nel",map.hdr.nel)
          write_dataset(h1,"noutflow",map.hdr.noutflow)
          write_dataset(h1,"npts",map.hdr.npts)
          write_dataset(h1,"nrank",map.hdr.nrank)
          write_dataset(h1,"ifsorted",ifsorted)

          # Data
          write_dataset(h2,"pmap",map.pmap)
          write_dataset(h2,"vmap",map.vmap)

          close(fid)
        end       # prank == nid0  

        return nothing
      end         # gen_rema2_h5
#---------------------------------------------------------------------- 
"""
      distribute_mesh(case::String, comm::MPI.Comm,nid0::Int)

      Read the case.re2 and case.ma2 files  and generate a case.rema2
      file with elements sorted according to the numbering in case.ma2.
      Boundary and Curve data also made consistent.

"""
      function distribute_mesh(casename::String,comm::MPI.Comm,nid0::Int)

        prank        = MPI.Comm_rank(comm)
        np          = MPI.Comm_size(comm)

        h5name      = casename*".rema2.h5"
        gnel        = 0
        ndim        = 0
        wdsize      = 4
        T           = Float32

        if prank == nid0
          fid      = h5open(h5name, "r")

          g        = fid["Re2"]
          g1       = g["Params"]
          g2       = g["Data"]

          h        = fid["Ma2"]
          h1       = h["Params"]
          h2       = h["Data"]

          gnel     = read(g1,"nelgt")
          ndim     = read(g1,"ndim")
          wdsize   = read(g1,"wdsize")

          # Read from the case.rema2.h5 file        
          xcg     = read(g2,"xc")
          ycg     = read(g2,"yc")
          if ndim == 3
            zcg   = read(g2,"zc")
          else
            zcg   = zeros(T,0,0)
          end

          pmapg   = read(h2,"pmap")
          vmapg   = read(h2,"vmap")

          close(fid)

        end  
      
        # Broadcast some header variables
        gnel       = MPI.bcast(gnel,      nid0, comm)
        ndim       = MPI.bcast(ndim,      nid0, comm)
        wdsize     = MPI.bcast(wdsize,    nid0, comm)

        if wdsize == 8
          T = Float64
        end  

        # Remaining Elements
        rem   = mod(gnel,np)  
        nel2  = gnel - rem
        lnel0 = floor(Int,gnel/np)
        lnel1 = lnel0+1
        # Add remaining elements to the end pranks
        last_prank = np-1     
        pextras = np - rem 

        i0    = floor(Int,log2(gnel/np))
        bsize = 2^i0

        # Allocate memory for xc,yc,pmap,vmap
        nc              = 2^ndim
        tupv            = fill(2,ndim+1)
        tupv[ndim+1]    = lnel1
        tup             = ntuple(i -> tupv[i],ndim+1)
        xc              = zeros(T,tup)
        yc              = zeros(T,tup)
        if ndim == 3
          zc  = zeros(T,tup)
        else
          zc  = zeros(T,0,0)
        end
       
        pmap  = fill(-1,lnel1)
        vmap  = fill(-1,tup)

        # Temporary array to store data for other processors        
        xcb    = zeros(T,nc,lnel1)
        ycb    = zeros(T,nc,lnel1)
        if ndim == 3
          zcb  = zeros(T,nc,lnel1)
        else
          zcb  = zeros(T,0,0)
        end
       
        pmapb  = fill(-1,lnel1)
        vmapb  = fill(-1,(nc,lnel1))

        if prank == nid0
          nid   = 0
          j     = 0
          elno  = pmapg[1] 
          while nid < np
            elst  = nid*bsize
            elend = (nid+1)*bsize
            i = 0
            while (elno>=elst && elno<elend) && j<gnel
              i     = i + 1
              j     = j + 1
              elno  = pmapg[j] 
              # Re2 Data              
              ind1  = CartesianIndices((1:nc,i:i))
              ind2  = CartesianIndices((1:nc,j:j)) 
              copyto!(xcb,ind1,xcg,ind2)
              copyto!(ycb,ind1,ycg,ind2)
              if ndim == 3
                copyto!(zcb,ind1,zcg,ind2)
              end  

              # Map data              
              copyto!(vmapb,ind1,vmapg,ind2)
              pmapb[i] = pmapg[j]
            end     # while j < jend
            nels = i
            # println("$nels $gnel $i $j $elno $elend")

            destn = nid
            if destn == nid0
              lnel  = nels
              # ind5  = CartesianIndices((1:nc,1:nels))
              npts  = nc*nels
              copyto!(xc,1,xcb,1,npts)
              copyto!(yc,1,ycb,1,npts)
              if ndim == 3
                copyto!(zc,1,zcb,1,npts)
              end  

              copyto!(vmap,1,vmapb,1,npts)
              # ind6  = CartesianIndices((1:nels))
              copyto!(pmap,1,pmapb,1,lnel)

            else  

              # Send no of elements
              send_status1 = MPI.Isend(nels,comm;dest=destn,tag=1)
             
              # Send Re2 data 
              send_status2 = MPI.Isend(xcb,comm;dest=destn,tag=2)
              send_status3 = MPI.Isend(ycb,comm;dest=destn,tag=3)
              if ndim == 3
                send_status4 = MPI.Isend(zcb,comm;dest=destn,tag=4)
              end    

              # Send ma2 data    
              send_status5 = MPI.Isend(vmapb,comm;dest=destn,tag=5)
              send_status6 = MPI.Isend(pmapb,comm;dest=destn,tag=6)

            end     # if nid == nid0
            nid = nid + 1
          end         # while nid < np 

        else  # prank != nid0

          nels = Vector{Int}(undef,1)
          # Receive nels
          recv_status1 = MPI.Irecv!(nels,comm; source=nid0,tag=1)

          # Receive Re2 data
          recv_status2 = MPI.Irecv!(xcb,comm; source=nid0,tag=2)
          recv_status3 = MPI.Irecv!(ycb,comm; source=nid0,tag=3)
          if ndim == 3
            recv_status4 = MPI.Irecv!(zcb,comm; source=nid0,tag=4)
          end  

          # Receive Ma2 data
          recv_status5 = MPI.Irecv!(vmapb,comm; source=nid0,tag=5)
          recv_status6 = MPI.Irecv!(pmapb,comm; source=nid0,tag=6)
          
          # Wait for data
          stat1 =  MPI.Wait!(recv_status1)
          stat2 =  MPI.Wait!(recv_status2)
          stat3 =  MPI.Wait!(recv_status3)
          if ndim == 3
            stat4 =  MPI.Wait!(recv_status4)
          end  
          stat5 =  MPI.Wait!(recv_status5)
          stat6 =  MPI.Wait!(recv_status6)

          # Move to Main arrays
          lnel  = nels[1]
          npts  = nc*lnel
          copyto!(xc,1,xcb,1,npts)
          copyto!(yc,1,ycb,1,npts)
          if ndim == 3
            copyto!(zc,1,zcb,1,npts)
          end  

          copyto!(vmap,1,vmapb,1,npts)
          # ind6  = CartesianIndices((1:nels))
          copyto!(pmap,1,pmapb,1,lnel)

        end   # if prank == nid0 

        MPI.Barrier(comm)

        tupv            = fill(2,ndim+1)
        tupv[ndim+1]    = lnel
        tup             = ntuple(i -> tupv[i],ndim+1)
        ind7            = CartesianIndices(tup)
        ind8            = CartesianIndices((1:lnel))

        # # Checks!        
        # # This should be commented out after checks        
        # # Write out an hdf5 file
        # ftest = casename*"_p$prank.h5"
        # fid2  = h5open(ftest, "w")
        # gg    = create_group(fid2,"Coords")
        # write_dataset(gg,"nel",lnel)

       
        # write_dataset(gg,"xc",xc[ind7])
        # write_dataset(gg,"yc",yc[ind7])
        # if ndim == 3
        #   write_dataset(gg,"zc",zc[ind7])
        # end
        # hh    = create_group(fid2,"Map")
        # write_dataset(hh,"pmap",pmap[ind8])
        # write_dataset(hh,"vmap",vmap[ind7])

        # close(fid2)
        # println("\n Writing Nel = $lnel, Elements, on Rank=$prank to $ftest\n")

        if ndim == 2
          return xc[ind7],yc[ind7],zc,pmap[ind8],vmap[ind7]
        else
          return xc[ind7],yc[ind7],zc[ind7],pmap[ind8],vmap[ind7]
        end
      end  # distribute_mesh 
#----------------------------------------------------------------------  
"""
      distribute_unsorted_mesh(case::String)

      Read the case.re2 and case.ma2 files and generate a case.rema2
      file. No Sorting of elements is done 

"""
      function distribute_unsorted_mesh(casename::String,comm::MPI.Comm,nid0::Int)

        prank        = MPI.Comm_rank(comm)
        np          = MPI.Comm_size(comm)

        h5name      = casename*".rema2.h5"
        gnel        = 0
        ndim        = 0
        wdsize      = 4
        T           = Float32

        if prank == nid0
          fid      = h5open(h5name, "r")

          g        = fid["Re2"]
          g1       = g["Params"]
          g2       = g["Data"]

          h        = fid["Ma2"]
          h1       = h["Params"]
          h2       = h["Data"]

          gnel     = read(g1,"nelgt")
          ndim     = read(g1,"ndim")
          wdsize   = read(g1,"wdsize")

          # Read from the case.rema2.h5 file        
          xcg     = read(g2,"xc")
          ycg     = read(g2,"yc")
          if ndim == 3
            zcg   = read(g2,"zc")
          else
            zcg   = zeros(T,0,0)
          end
 
          pmapg   = read(h2,"pmap")
          vmapg   = read(h2,"vmap")

          key     = sortperm(pmapg) # Sorting index

          close(fid)

        end  

        # Broadcast some header variables        
        gnel       = MPI.bcast(gnel,      nid0, comm)
        ndim       = MPI.bcast(ndim,      nid0, comm)
        wdsize     = MPI.bcast(wdsize,    nid0, comm)

        if wdsize == 8
          T = Float64
        end  

        # Remaining Elements
        rem   = mod(gnel,np)  
        nel2  = gnel - rem
        lnel0 = floor(Int,nel2/np)
        lnel1 = lnel0+1
        # Add remaining elements to the end pranks
        last_prank = np-1     
        pextras = np - rem 

        i0    = floor(Int,log2(gnel/np))
        bsize = 2^i0

        # Allocate memory for xc,yc,pmap,vmap
        nc    = 2^ndim
        xc    = zeros(T,nc,lnel1)
        yc    = zeros(T,nc,lnel1)
        if ndim == 3
          zc  = zeros(T,nc,lnel1)
        else
          zc  = zeros(T,0,0)
        end
        pmap  = fill(-1,lnel1)
        vmap  = fill(-1,(nc,lnel1))

        # Temporary array to store data for other processors        
        xcb    = zeros(T,nc,lnel1)
        ycb    = zeros(T,nc,lnel1)
        if ndim == 3
          zcb  = zeros(T,nc,lnel1)
        else
          zcb  = zeros(T,0,0)
        end
       
        pmapb  = fill(-1,lnel1)
        vmapb  = fill(-1,(nc,lnel1))

        if prank == nid0
          nid   = 0
          j     = 0
          k     = key[1]
          elno  = pmapg[k] 
          while nid < np

            elst  = nid*bsize
            elend = (nid+1)*bsize
            i = 0
            while (elno>=elst && elno<elend) && j<gnel
              i     = i + 1
              j     = j + 1
              k     = key[j]
              elno  = pmapg[k] 
             
              # Re2 Data              
              ind1  = CartesianIndices((1:nc,i:i))
              ind2  = CartesianIndices((1:nc,k:k)) 
              copyto!(xcb,ind1,xcg,ind2)
              copyto!(ycb,ind1,ycg,ind2)
              if ndim == 3
                copyto!(zcb,ind1,zcg,ind2)
              end  

              # Map data              
              copyto!(vmapb,ind1,vmapg,ind2)
              pmapb[i] = pmapg[k]

            end     # while j < jend
            nels = i

            destn  = nid
            if destn == nid0
              lnel  = nels
              ind5  = CartesianIndices((1:nc,1:nels))
              copyto!(xc,ind5,xcb,ind5)
              copyto!(yc,ind5,ycb,ind5)
              if ndim == 3
                copyto!(zc,ind1,zcb,ind2)
              end  

              copyto!(vmap,ind5,vmapb,ind5)
              ind6  = CartesianIndices((1:nels))
              copyto!(pmap,ind6,pmapb,ind6)

            else  

              # Send no of elements
              send_status1 = MPI.Isend(nels,comm;dest=destn,tag=1)
             
              # Send Re2 data 
              send_status2 = MPI.Isend(xcb,comm;dest=destn,tag=2)
              send_status3 = MPI.Isend(ycb,comm;dest=destn,tag=3)
              if ndim == 3
                send_status4 = MPI.Isend(zcb,comm;dest=destn,tag=4)
              end    

              # Send ma2 data    
              send_status5 = MPI.Isend(vmapb,comm;dest=destn,tag=5)
              send_status6 = MPI.Isend(pmapb,comm;dest=destn,tag=6)

            end     # if nid == nid0
            nid = nid + 1
          end         # while nid < np 

        else  # prank != nid0

          nels = Vector{Int}(undef,1)
          # Receive nels
          recv_status1 = MPI.Irecv!(nels,comm; source=nid0,tag=1)

          # Receive Re2 data
          recv_status2 = MPI.Irecv!(xcb,comm; source=nid0,tag=2)
          recv_status3 = MPI.Irecv!(ycb,comm; source=nid0,tag=3)
          if ndim == 3
            recv_status4 = MPI.Irecv!(zcb,comm; source=nid0,tag=4)
          end  

          # Receive Ma2 data
          recv_status5 = MPI.Irecv!(vmapb,comm; source=nid0,tag=5)
          recv_status6 = MPI.Irecv!(pmapb,comm; source=nid0,tag=6)
          
          # Wait for data
          stat1 =  MPI.Wait!(recv_status1)
          stat2 =  MPI.Wait!(recv_status2)
          stat3 =  MPI.Wait!(recv_status3)
          if ndim == 3
            stat4 =  MPI.Wait!(recv_status4)
          end  
          stat5 =  MPI.Wait!(recv_status5)
          stat6 =  MPI.Wait!(recv_status6)

          # Move to Main arrays
          lnel  = nels[1]
          ind1  = CartesianIndices((1:nc,1:lnel))
          copyto!(xc,ind1,xcb,ind1)
          copyto!(yc,ind1,ycb,ind1)
          if ndim == 3
            copyto!(zc,ind1,zcb,ind2)
          end  

          copyto!(vmap,ind1,vmapb,ind1)
          ind2  = CartesianIndices((1:lnel))
          copyto!(pmap,ind2,pmapb,ind2)

        end   # if prank == nid0 

        MPI.Barrier(comm)

        # Checks!        
        # This should be commented out after checks        
        # Write out an hdf5 file
        ftest = casename*"_p$prank.h5"
        fid2  = h5open(ftest, "w")
        gg    = create_group(fid2,"Coords")
        write_dataset(gg,"nel",lnel)
        ind7  = CartesianIndices((1:nc,1:lnel))
        write_dataset(gg,"xc",xc[ind7])
        write_dataset(gg,"yc",yc[ind7])
        if ndim == 3
          write_dataset(gg,"zc",zc[ind7])
        end
        hh    = create_group(fid2,"Map")
        ind8  = CartesianIndices((1:lnel))
        write_dataset(hh,"pmap",pmap[ind8])
        write_dataset(hh,"vmap",vmap[ind7])

        close(fid2)

        println("\n Writing Nel = $lnel, Elements, on Rank=$prank to $ftest\n")

        if ndim == 2
          return xc[ind7],yc[ind7],zc,pmap[ind8],vmap[ind7]
        else
          return xc[ind7],yc[ind7],zc[ind7],pmap[ind8],vmap[ind7]
        end
        # return xc,yc,zc,pmap,vmap
      end  # distribute_unsorted_mesh 
#----------------------------------------------------------------------  

"""
  Convert from Preprocessor Coordinates indices 
  to Symmetric Coordinate indicies

  Preprocessor Corner notation  ------->  Symmetric Corner notation:

          4+-----+3    ^ s                    3+-----+4    ^ s
          /     /|     |                      /     /|     |
         /     / |     |                     /     / |     |
       8+-----+7 +2    +----> r            7+-----+8 +2    +----> r
        |     | /     /                     |     | /     /
        |     |/     /                      |     |/     /
       5+-----+6    t                      5+-----+6    t


  preproctosymm(i::Int)

# Examples
```julia-repl
julia> j = preproctosymm(3)
4
```
"""
function preproctosymm(i::Int)

#  indx = [1; 2; 4; 3; 5; 6; 8; 7]

  si = i
  if i==3
    si = 4
  elseif i==4
    si = 3
  elseif i==7
    si = 8
  elseif i==8
    si = 7
  end

  if i<1 || i>8
    error("preproctosymm: Index Out of bounds: $i")
  end

  return si
end  

#----------------------------------------------------------------------

"""
  Convert from Preprocessor Coordinates indices 
  to Symmetric Coordinate indicies

  Symmetric Corner notation: ------->  Preprocessor Corner notation    
                             
         3+-----+4    ^ s                     4+-----+3    ^ s                
         /     /|     |                       /     /|     |                  
        /     / |     |                      /     / |     |                  
      7+-----+8 +2    +----> r             8+-----+7 +2    +----> r           
       |     | /     /                      |     | /     /                   
       |     |/     /                       |     |/     /                    
      5+-----+6    t                       5+-----+6    t                     


  pymmtopreproc(i::Int)

# Examples
```julia-repl
julia> j = symmtopreproc(7)
8
```
"""
function symmtopreproc(i::Int)

#  indx = [1; 2; 4; 3; 5; 6; 8; 7]
  
  ppi = i
  if i==3
    ppi = 4
  elseif i==4
    ppi = 3
  elseif i==7
   ppi = 8
  elseif i==8
   ppi = 7
  end

  if i<1 || i>8
    error("symmtopreproc: Index Out of bounds: $i")
  end

  return ppi
end  
#----------------------------------------------------------------------
@doc raw"""

      edgeindices(ei::Int,dims::Int)

      Return Tuple of vertex numbers for edge index ei.

# Examples
```julia-repl
julia> edge2d = edgeindices(3,2)
((1,1), (1,2))
julia> edge3d = edgeindices(7,3)
((1,1,2), (1,2,2))
```

"""
function edgeindices(ei::Int,dims::Int)

    if dims == 1
      tup = tuple()
    elseif dims == 2
      tup = edgeindices2d(ei)
    elseif dims == 3
      tup = edgeindices3d(ei)
    else
      error("edge indices for dim>3 not defined: $dims")
    end

  return tup
end
#---------------------------------------------------------------------- 

@doc raw"""

      edgeindices2d

      Symmetric Edge numbering {i} in 
      Symmetric Corner notation

      (1,2)                (2,2)
       3+--------{2}--------+4       ^ s       
        |                   |        |
        |                   |        |  
        |                   |        |  
        |                   |        |         
       {3}                 {4}       |        
        |                   |        | 
        |                   |        |
        |                   |        |
        |                   |        |
       1+--------{1}--------+2       +--------------> r  
      (1,1)               (2,1)


"""
function edgeindices2d(ei::Int)

  @assert ei>0 && ei<=4 "Unknown edge index: $ei"

  if ei == 1
    t1    = tuple(1,1)
    t2    = tuple(2,1)
  elseif ei == 2
    t1    = tuple(1,2)
    t2    = tuple(2,2)
  elseif ei == 3
    t1    = tuple(1,1)
    t2    = tuple(1,2)
  elseif ei == 4 
    t1    = tuple(2,1)
    t2    = tuple(2,2)
  end

  return t1,t2
end
#---------------------------------------------------------------------- 
@doc raw"""

      edgeindicese3d

      Symmetric Edge numbering {i} in 
      Symmetric Corner notation

                   (1,2,2)                  (2,2,2)
                      7+---------{4}----------+8                ^ s       
                      /|                     /|                 |  
                     /                      / |                 | 
                  {11} |                   /  |                 |
                   /  {7}                {12} |                 |
                  /    |                 /    |                 |  
        (1,2,1) 3+--------{2}-----------+4   {8}                |         
                 |     |                |     |                 |        
                 |                      |     |                 |
                 |     |(1,1,2)         |     |                 |
                 |    5+- - - - {3} - - |- - -+6 (2,1,2)        +---------------> r
                {5}   /                {6}   /                 /
                 |   /                  |   /                 /
                 | {9}                  | {10}               /
                 | /                    | /                 / 
                 |/                     |/                 / 
                1+---------{1}----------+2                /  
             (1,1,1)                 (2,1,1)             t


"""
function edgeindices3d(ei::Int)

  @assert ei>0 && ei<=12 "Unknown edge index: $ei"

  if ei == 1
    t1    = tuple(1,1,1)
    t2    = tuple(2,1,1)
  elseif ei == 2
    t1    = tuple(1,2,1)
    t2    = tuple(2,2,1)
  elseif ei == 3
    t1    = tuple(1,1,2)
    t2    = tuple(2,1,2)
  elseif ei == 4 
    t1    = tuple(1,2,2)
    t2    = tuple(2,2,2)
  elseif ei == 5
    t1    = tuple(1,1,1)
    t2    = tuple(1,2,1)
  elseif ei == 6
    t1    = tuple(2,1,1)
    t2    = tuple(2,2,1)
  elseif ei == 7
    t1    = tuple(1,1,2)
    t2    = tuple(1,2,2)
  elseif ei == 8
    t1    = tuple(2,1,2)
    t2    = tuple(2,2,2)
  elseif ei == 9
    t1    = tuple(1,1,1)
    t2    = tuple(1,1,2)
  elseif ei == 10
    t1    = tuple(2,1,1)
    t2    = tuple(2,1,2)
  elseif ei == 11
    t1    = tuple(1,2,1)
    t2    = tuple(1,2,2)
  elseif ei == 12
    t1    = tuple(2,2,1)
    t2    = tuple(2,2,2)
  end

  return t1,t2
end
#---------------------------------------------------------------------- 






# testing the re2/map/par readers
#

#      using Revise
#      using PyCall
#      using PyPlot 
      using HDF5
      
      include("SEMla.jl")            # SEMla
#      include("NekTools.jl")
#      include("$(JULIACOMMON)/MoveFigure.jl")
      
      using MPI

      if !MPI.Initialized()
        MPI.Init()
        const comm      = MPI.COMM_WORLD
        const rank      = MPI.Comm_rank(comm)
        const comm_size = MPI.Comm_size(comm)
      end  

#     Node at which we want to read rema2 file      
      nid0  = 0

      casename = "cyl"

#     For now we just make sure we generate the .rema2.h5 file.      
      SEMla.gen_rema2(casename,nid0,comm)

      h5name      = casename*".rema2.h5"
      gnel        = 0
      ndim        = 0
      wdsize      = 4
      T           = Float32
      
      if rank == nid0
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

#        xc       = read(g2,"xc")
#        yc       = read(g2,"yc")

#        pmap     = read(h2,"pmap")    # Processor Map
#        vmap     = read(h2,"vmap")    # Vertex Map
      end  

      gnel        = MPI.bcast(gnel, nid0, comm)
      ndim        = MPI.bcast(ndim, nid0, comm)
      wdsize      = MPI.bcast(wdsize, nid0, comm)

      if wdsize == 8
        T = Float64
      end  

      println("\n Nel = $gnel, Ndim=$ndim, wdsize=$wdsize, on Rank=$rank\n")

      np    = comm_size    # No of Processors
      rem   = mod(gnel,np)
      nel2  = gnel - rem
      lnel0 = floor(Int,nel2/np)
      # Add remaining elements to the end ranks
      last_rank = np-1     
      pextras = np - rem 
   
      lnel  = lnel0
      if rank >= pextras
        lnel = lnel0 + 1
      end

      # Allocate memory for xc,yc
      nc    = 2^ndim
      xc    = zeros(T,nc,lnel)
      yc    = zeros(T,nc,lnel)

      # Temporary array to store data for other processors        
      xcb   = zeros(T,nc,lnel0+1)
      ycb   = zeros(T,nc,lnel0+1)

      if rank == nid0
        # Read from the case.rema2.h5 file        
        xcg  = read(g2,"xc")
        ycg  = read(g2,"yc")
       
        j    = 0
        nid  = 0
        jend = 0
        while nid < np
          global jend, j, nid
          jend  = jend + lnel0
          if nid >= pextras
            jend   = jend + 1
          end

          i    = 0
          while j < jend
            j = j + 1
            i = i + 1

            ind1  = CartesianIndices((1:nc,i:i))
            ind2  = CartesianIndices((1:nc,j:j)) 
            copyto!(xcb,ind1,xcg,ind2)
            copyto!(ycb,ind1,ycg,ind2)

          end     # while j < jend

          if nid == nid0
            ind1  = CartesianIndices((1:nc,1:lnel))
            copyto!(xc,ind1,xcb,ind1)
            copyto!(yc,ind1,ycb,ind1)
          else  
#            vwx          = view(xcin[:,:,1:jend]) 
            destn        = nid
            send_status1 = MPI.Isend(xcb,comm;dest=destn)

#            vwy          = view(ycb[:,:,1:jend]) 
            destn        = nid
            send_status2 = MPI.Isend(ycb,comm;dest=destn)
          end     # if nid == nid0

          nid = nid + 1

        end         # while nid < np 

      else  # rank != nid0
         recv_status1 = MPI.Irecv!(xcb,comm; source=nid0)
         recv_status2 = MPI.Irecv!(ycb,comm; source=nid0)

         stat1 =  MPI.Wait!(recv_status1)
         stat2 =  MPI.Wait!(recv_status2)

         ind3  = CartesianIndices((1:nc,1:lnel))
         copyto!(xc,ind3,xcb,ind3)
         copyto!(yc,ind3,ycb,ind3)
      end   # if rank == nid0 

      MPI.Barrier(comm)

#     Write out an hdf5 file
      ftest = casename*"_p$rank.h5"
      fid2  = h5open(ftest, "w")
      gg    = create_group(fid2,"Coords")
      write_dataset(gg,"xc",xc)
      write_dataset(gg,"yc",yc)
      close(fid2)

      println("\n Writing Nel = $lnel, Elements, on Rank=$rank to $ftest\n")

#     Distribute Mesh across Processors
#--------------------------------------------------     

            




#      length   = nc*nel
#      send_buf = rand(Float64,length) 
#      recv_buf = Vector{Float64}(undef,length)
#
#      destn = (rank+1)%wsize
#      send_status = MPI.Isend(send_buf,comm;dest=destn)
#
#      src = rank-1
#      if src<0
#        src = wsize-1
#      end
#      recv_status = MPI.Irecv!(recv_buf,comm; source=src)
#
#      stat =  MPI.Wait!(recv_status)

#-------------------------------------------------- 

      if rank == nid0
        close(fid)
        println("Done")
      end  

      MPI.Finalize()

#      if rank == 0
#        println("Done. MPI not Finalized")
#      end  


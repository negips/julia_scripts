# testing the re2/map/par readers
#

#      using Revise
#      using PyCall
#      using PyPlot 
      using HDF5
      
      include("SEMla.jl")            # SEMla
      include("NekTools.jl")
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
      nel         = 0
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

      println("\n Nel = $nelg, Ndim=$ndim, wdsize=$wdsize, on Rank=$rank\n")

      np    = comm_size    # No of Processors
      rem   = mod(gnel,np)
      nel2  = gnel - rem
      lnel0 = floor(Int,nel2/np)
#     Add remaining elements to the end ranks
      last_rank = np-1     
      pextras = np - rem 
   
      lnel  = lnel0
      if rank >= pextras
        lnel = lnel0 + 1
      end
      xc    = zeros(T,nc,lnel) 
      yc    = zeros(T,nc,lnel) 

#     Allocate memory for xc,yc
      nc    = 2^ndim
      xc    = zeros(T,nc,lnel)
      yc    = zeros(T,nc,lnel)

      if rank != nid0
#         recv_status1 = MPI.Irecv!(xc,comm; source=nid0)
#         recv_status2 = MPI.Irecv!(xc,comm; source=nid0)
      else
        xcg      = read(g2,"xc")
        ycg      = read(g2,"yc")
        
        for i in 0:last_rank
          if rank == nid
          end
        end
      end  

      MPI.Barrier(comm)


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


# testing the re2/map/par readers
#

#      using Revise
#      using PyCall
#      using PyPlot 
      using HDF5
      
      include("JNek_IO.jl")            # JNek_IO
      include("NekTools.jl")
#      include("$(JULIACOMMON)/MoveFigure.jl")
      
      using MPI

      if !MPI.Initialized()
        MPI.Init()
        const comm = MPI.COMM_WORLD
        const rank = MPI.Comm_rank(comm)
      end  

#     Node at which we want to read rema2 file      
      nid0  = 0

      casename = "cyl"

#     For now we just make sure we generate the .rema2.h5 file.      
      JNek_IO.gen_rema2(casename,nid0,comm)

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

        nel      = read(g1,"nelgt")
        ndim     = read(g1,"ndim")
        wdsize   = read(g1,"wdsize")

        xc       = read(g2,"xc")
        yc       = read(g2,"yc")

        pmap     = read(h2,"pmap")    # Processor Map
        vmap     = read(h2,"vmap")    # Vertex Map
      end  

      nel         = MPI.bcast(nel,  nid0, comm)
      ndim        = MPI.bcast(ndim, nid0, comm)
      wdsize      = MPI.bcast(wdsize, nid0, comm)

      if wdsize == 8
        T = Float64
      end  

      println("\n Nel = $nel, Ndim=$ndim, wdsize=$wdsize, on Rank=$rank\n")

#     Allocate memory for xc,yc
      nc    = 2^ndim
#      xc    = zeros()

      MPI.Barrier(comm)

      if rank == nid0
        close(fid)
        println("Done")
      end  

      MPI.Finalize()

#      if rank == 0
#        println("Done. MPI not Finalized")
#      end  


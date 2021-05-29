# testing the re2/map/par readers
#


#      include("JNek_PARALLEL.jl")         # JNek_PARALLEL
      include("JNek_IO.jl")               # JNek_IO

      using .JNek_IO
      
      using MPI

      nid0  = 0

      if (MPI.Initialized() == false)
        MPI.Init()
      end  
        
      comm = MPI.COMM_WORLD
      rank = MPI.Comm_rank(comm)

      re2   = "channelp.re2"
      map   = "channelp.ma2"

      hdr,version,nelgt,ldimr,nelgv = read_re2(re2,nid0)

      if rank == 0
        println("Done")
      end  

      MPI.Finalize()

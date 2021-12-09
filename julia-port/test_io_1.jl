# testing the re2/map/par readers
#


#      include("JNek_MPI.jl")         # JNek_MPI
      include("JNek_IO_1.jl")          # JNek_IO
      
      using MPI

      if (MPI.Initialized() == false)
        println("MPI = $(MPI.Initialized())")
        MPI.Init()
      end  

      re2   = "channelp.re2"
      map   = "channelp.ma2"
      fld   = "taylor0.f00001"

      nid0  = 0
      comm  = MPI.COMM_WORLD

      if MPI.Comm_rank(comm) == nid0
        println("Started MPI with $(MPI.Comm_size(comm)) ranks\n")
      end

      hdr,version,nelgt,ldimr,nelgv = JNek_IO.read_re2(re2,MPI,nid0)
      
      hdr,version,wdsize,nx,ny,nz,nel,nelgt,time,istep,fid0,nfileo,rdcode,p0th,ifprmesh = JNek_IO.read_fld(fld,MPI,nid0)

      println("Done")

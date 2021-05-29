#     Port for reader_par.f
#     Author:     Prabal Negi
#

      module JNek_IO

      using MPI

      export read_re2_hdr,
             read_re2, 
             read_ma2

      function __init__()

        if MPI.Initialized() == false

          MPI.Init()

        end    
          
        comm = MPI.COMM_WORLD
        rank = MPI.Comm_rank(comm)

        if rank == 0
          println("Initialied MPI in Module JNek_IO")
        end  

        return nothing
      end 

#----------------------------------------------------------------------

      function read_re2_hdr(fid::IOStream, rank)

#        hdr = repeat(" ", 26)
#        version = repeat(" ", 5)
#        nelgt = 0
#        nelgv = 0
#        ldimr = 0

        println("Reading Header on rank $(rank)")
        
        nbytes  = 26
        hdrutf  = read(fid,nbytes)
        hdrC    = Char.(hdrutf)
        version = String(hdrC[1:5])
        nelgtS  = String(hdrC[6:14])
        ldimrS  = String(hdrC[15:17])
        nelgvS  = String(hdrC[18:26])

        nelgt   = parse(Int, nelgtS)
        ldimr   = parse(Int, ldimrS)
        nelgv   = parse(Int, nelgvS)

        hdr     = String(hdrC) 


        return hdr,version,nelgt,ldimr,nelgv
      end     # read_re2_hdr

#---------------------------------------------------------------------- 

      function read_re2(f::String, nid0::Int64)

        hdr = repeat(" ", 26)
        version = repeat(" ", 5)
        nelgt = 0
        nelgv = 0
        ldimr = 0

        comm = MPI.COMM_WORLD

        if MPI.Comm_rank(comm) == nid0
          println("Reading $(f) on rank $(MPI.Comm_rank(comm))\n")
        end
        
        MPI.Barrier(comm)
        rank = MPI.Comm_rank(comm)

        if (rank == nid0)
          fid = open(f, "r")
        end  

        if rank == nid0
          hdr,version,nelgt,ldimr,nelgv = read_re2_hdr(fid,rank)
        end  


        buf = MPI.Buffer(hdr,26,MPI.CHAR)
        MPI.Bcast!(buf,     nid0,comm)

        buf = MPI.Buffer(version,5,MPI.CHAR)
        MPI.Bcast!(buf,     nid0,comm)

        nelgt = MPI.bcast(nelgt,     nid0,comm)
        ldimr = MPI.bcast(ldimr,     nid0,comm)
        nelgv = MPI.bcast(nelgv,     nid0,comm)

#       Read the data here        

        if rank == nid0
          close(fid)
        end  

        return hdr,version,nelgt,ldimr,nelgv
      end     # read_re2

#---------------------------------------------------------------------- 

      function read_ma2(f::String, nid0::Int64)


        comm = MPI.COMM_WORLD

        if MPI.Comm_rank(comm) == nid0
          println("Reading $(f) on rank $(MPI.Comm_rank(comm))")
        end
        
        MPI.Barrier(comm)

        rank = MPI.Comm_rank(comm)

        if rank == nid0
          fid = open(f, "r")
        end

        close(fid)

        return nothing
      end     # read_ma2

#---------------------------------------------------------------------- 

      end   # Module JNek_IO












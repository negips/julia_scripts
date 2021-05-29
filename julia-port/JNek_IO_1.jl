#     Port for reader_par.f
#     Author:     Prabal Negi
#

      module JNek_IO

      export read_re2_hdr,
             read_re2, 
             read_ma2

      function read_re2_hdr(fid::IOStream, rank, nid0)

#        comm = MPI.COMM_WORLD
        
        hdr     = repeat(" ", 26)
        version = repeat(" ", 5)
        nelgt   = 0
        nelgv   = 0
        ldimr   = 0

        if rank == nid0
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
        end


        return hdr,version,nelgt,ldimr,nelgv
      end     # read_re2_hdr

#---------------------------------------------------------------------- 

      function read_re2(f::String, MPI::Module, nid0::Int64)


        comm = MPI.COMM_WORLD

        if MPI.Comm_rank(comm) == nid0
          println("Reading $(f) on rank $(MPI.Comm_rank(comm))\n")
        end
        
        MPI.Barrier(comm)

        fid = open(f, "r")

        rank = MPI.Comm_rank(comm)

        hdr,version,nelgt,ldimr,nelgv = read_re2_hdr(fid,rank, nid0);

        buf = MPI.Buffer(hdr,length(hdr),MPI.CHAR)
        MPI.Bcast!(buf,     nid0,comm)

        buf = MPI.Buffer(version,length(version),MPI.CHAR)
        MPI.Bcast!(buf,     nid0,comm)

        nelgt = MPI.bcast(nelgt,     nid0,comm)
        ldimr = MPI.bcast(ldimr,     nid0,comm)
        nelgv = MPI.bcast(nelgv,     nid0,comm)

#       Read the data here        

        close(fid)

        return hdr,version,nelgt,ldimr,nelgv
      end     # read_re2

#---------------------------------------------------------------------- 

      function read_ma2(f::String, MPI::Module, nid0::Int64)


        comm = MPI.COMM_WORLD

        if MPI.Comm_rank(comm) == nid0
          println("Reading $(f) on rank $(MPI.Comm_rank(comm))")
        end
        
        MPI.Barrier(comm)

        fid = open(f, "r")

        rank = MPI.Comm_rank(comm)

#        MPI.Bcast!(hdr,nid0,MPI.CHAR)
#        MPI.Bcast!(version,   nid0,MPI.CHAR)
#        MPI.Bcast!(nelgt,     nid0,MPI.INT64_T)
#        MPI.Bcast!(ldimr,     nid0,MPI.INT64_T)
#        MPI.Bcast!(ldimr,     nid0,MPI.INT64_T)

#       testing        
        if rank == 1
          println(hdr)
        end  

        close(fid)

#        MPI.Finalize()

        return hdr,version,nelgt,ldimr,nelgv
      end     # read_ma2

#---------------------------------------------------------------------- 

      end   # Module JNek_IO












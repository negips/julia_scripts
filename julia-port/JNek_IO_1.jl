#     Port for reader_par.f
#     Author:     Prabal Negi
#

      module JNek_IO

      export read_re2_hdr,
             read_re2, 
             read_ma2,
             read_fld_std_hdr

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
      function read_fld(f::String, MPI::Module, nid0::Int64)


        comm = MPI.COMM_WORLD

        if MPI.Comm_rank(comm) == nid0
          println("Reading $(f) on rank $(MPI.Comm_rank(comm))\n")
        end
        
        MPI.Barrier(comm)

        fid = open(f, "r")

        rank = MPI.Comm_rank(comm)

        hdr,version,wdsize,nx,ny,nz,nel,nelgt,time,istep,fid0,nfileo,rdcode,p0th,ifprmesh = read_fld_std_hdr(fid,rank,nid0)

        buf = MPI.Buffer(hdr,length(hdr),MPI.CHAR)
        MPI.Bcast!(buf,     nid0,comm)

        buf = MPI.Buffer(version,length(version),MPI.CHAR)
        MPI.Bcast!(buf,     nid0,comm)

        wdsize          = MPI.bcast(wdsize,     nid0,comm)
        nx              = MPI.bcast(nx,         nid0,comm)
        ny              = MPI.bcast(ny,         nid0,comm)
        nz              = MPI.bcast(nz,         nid0,comm)
        nel             = MPI.bcast(nel,        nid0,comm)
        nelgt           = MPI.bcast(nelgt,      nid0,comm)
        time            = MPI.bcast(time,       nid0,comm)
        istep           = MPI.bcast(istep,      nid0,comm)
        fid0            = MPI.bcast(fid0,       nid0,comm)
        nfileo          = MPI.bcast(nfileo,     nid0,comm)

        buf             = MPI.Buffer(rdcode,length(rdcode),MPI.CHAR)
        MPI.Bcast!(buf,     nid0,comm)
        p0th            = MPI.bcast(p0th,       nid0,comm)
        ifprmesh        = MPI.bcast(ifprmesh,   nid0,comm)


#       Read the data here        

        close(fid)

        return hdr,version,wdsize,nx,ny,nz,nel,nelgt,time,istep,fid0,nfileo,rdcode,p0th,ifprmesh
      end     # read_fld

#---------------------------------------------------------------------- 

      function read_fld_std_hdr(fid::IOStream, rank, nid0)

#        comm = MPI.COMM_WORLD
        
        hdr           = repeat(" ", 132)
        version       = repeat(" ", 5)
        wdsize        = 0 
        nx            = 0
        ny            = 0
        nz            = 0
        nel           = 0
        nelgt         = 0
        time          = 0.0
        istep         = 0
        fid0          = 0
        nfileo        = 0
        rdcode        = repeat(" ",10)
        p0th          = 0.0 
        ifprmesh      = false

        if rank == nid0
          println("Reading Header on rank $(rank)")
        
          nbytes        = 132
          hdrutf        = read(fid,nbytes)
          hdrC          = Char.(hdrutf)
          
          st            = 2         # step
          i             = 1
          j             = 4
          version       = String(hdrC[i:j])
          i             = j+st # 7 
          j             = i+0
          wdsizS        = String(hdrC[i:j])
          i             = j+st # 
          j             = i+1
          nxS           = String(hdrC[i:j])
          i             = j+st # 
          j             = i+1
          nyS           = String(hdrC[i:j])
          i             = j+st # 
          j             = i+1
          nzS           = String(hdrC[i:j])
          i             = j+st # 
          j             = i+9
          nelS          = String(hdrC[i:j])
          i             = j+st # 
          j             = i+9
          nelgS         = String(hdrC[i:j])
          i             = j+st # 
          j             = i+19
          timeS         = String(hdrC[i:j])
          i             = j+st # 
          j             = i+8
          istepS        = String(hdrC[i:j])
          i             = j+st # 
          j             = i+5
          fid0S         = String(hdrC[i:j])
          i             = j+st # 
          j             = i+5
          nfileoS       = String(hdrC[i:j])
          i             = j+st # 
          j             = i+9
          rdcodeS       = String(hdrC[i:j])
          i             = j+st # 
          j             = i+14
          p0thS         = String(hdrC[i:j])
          i             = j+st # 
          j             = i+0
          ifpr_meshS    = String(hdrC[i:j])
    
          wdsize        = parse(Int, wdsizS)
          nx            = parse(Int, nxS)
          ny            = parse(Int, nyS)
          nz            = parse(Int, nzS)
          nel           = parse(Int, nelS)
          nelgt         = parse(Int, nelgS)
          time          = parse(Float64, timeS)
          istep         = parse(Int, istepS)
          fid0          = parse(Int, fid0S)
          nfileo        = parse(Int, nfileoS)
          rdcode        = rdcodeS
          p0th          = parse(Float64, p0thS)
          if (ifpr_meshS == "T")
            ifprmesh = true
          else
            ifprmesh = false
          end  
          hdr           = String(hdrC) 
        end


        return hdr,version,wdsize,nx,ny,nz,nel,nelgt,time,istep,fid0,nfileo,rdcode,p0th,ifprmesh
      end     # read_fld_std_hdr

#---------------------------------------------------------------------- 

      end   # Module JNek_IO












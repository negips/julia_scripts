#     Author:     Prabal Negi
#

      module JNek_IO

      include("JNek_IO_Abstract.jl")
      include("JNek_IO_Structs.jl")
      include("JNek_IO_Constructors.jl")
      include("JNek_IO_Extends.jl")

      using MPI
      using HDF5


      export NekField,
             Re2Field,
             TwoTensorField


      export read_re2_hdr,
             read_re2,
             read_re2_struct,
             read_ma2_hdr,
             read_ma2,
             read_fld,
             read_fld_struct,
             gen_rema2

#----------------------------------------------------------------------  

      function __init__()

#        if MPI.Initialized() == false
#
#          MPI.Init()
#
#        end    
#          
#        comm = MPI.COMM_WORLD
#        rank = MPI.Comm_rank(comm)
#
#        if rank == 0
#          println("Initialied MPI in Module JNek_IO")
#        end  

        return nothing
      end 

#----------------------------------------------------------------------

      function init()

#       Can't reinitialize MPI        

        return nothing
      end 

#----------------------------------------------------------------------

"""
      read_re2_hdr(fid::IOStream)

      Read the header of IOStream from a .re2 file.

"""
      function read_re2_hdr(fid::IOStream)

        println(".re2: Reading Header")

        nbytes  = 20*4
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

        test    = read(fid,Float32)

        if_byte_swap = byte_swap_test(test)

        return hdr,version,nelgt,ldimr,nelgv,if_byte_swap
      end     # read_re2_hdr

#----------------------------------------------------------------------
"""
      read_re2(f::String, nid0::Int64, comm::MPI.Comm)

      Read the re2 file (f) at MPI process ID=nid0 over the
      MPI Communicator Comm.

      Output: wdsizi,hdr,version,nelgt,ldimr,nelgv,xc,yc,zc,ncurve,
              curveieg,curveiside,curveparam,curvetype,cbl,bl

"""
      function read_re2(f::String, nid0::Int64, comm::MPI.Comm)

        hdr = repeat(" ", 26)
        version = repeat(" ", 5)
        nelgt = 0
        nelgv = 0
        ldimr = 0

        rank = MPI.Comm_rank(comm)

        if rank == nid0
          println("Reading $(f) on rank $rank")

          fid = open(f, "r")
          hdr,version,nelgt,ldimr,nelgv,if_byte_swap = read_re2_hdr(fid)
        end  

        MPI.Barrier(comm)

        buf = MPI.Buffer(hdr,26,MPI.CHAR)
        MPI.Bcast!(buf,     nid0,comm)

        buf = MPI.Buffer(version,5,MPI.CHAR)
        MPI.Bcast!(buf,     nid0,comm)

        nelgt = MPI.bcast(nelgt,     nid0,comm)
        ldimr = MPI.bcast(ldimr,     nid0,comm)
        nelgv = MPI.bcast(nelgv,     nid0,comm)

        wdsizi = 4
        if cmp(version,"#v002") == 0 || cmp(version, "#v003")
          wdsizi = 8
        end   

#       Read the Mesh data here
        xc,yc,zc = read_re2_mesh(fid,nid0,ldimr,nelgt,wdsizi,comm)

        ncurve,curveieg,curveiside,curveparam,curvetype = read_re2_curve(fid,nid0,ldimr,nelgt,wdsizi)

        cbl,bl = read_re2_bc(fid,nid0,ldimr,nelgt,wdsizi)

        if rank == nid0
          close(fid)
        end 

        return wdsizi,hdr,version,nelgt,ldimr,nelgv,xc,yc,zc,ncurve,curveieg,curveiside,curveparam,curvetype,cbl,bl
      end     # read_re2

#---------------------------------------------------------------------- 
"""
      read_re2_struct(f::String, nid0::Int64, comm::MPI.Comm)

      Read the re2 file (f) at MPI process ID=nid0 over the
      MPI Communicator Comm and output structured data Re2Field

      Output: Re2Field 

"""
      function read_re2_struct(f::String, nid0::Int64, comm::MPI.Comm)

        re2fld = Re2Field()

        rank = MPI.Comm_rank(comm)

        if rank == nid0
          println("Reading $(f) on rank $rank")

          fid = open(f, "r")
          hdr,version,nelgt,ldimr,nelgv,if_byte_swap = read_re2_hdr(fid)

          wdsizi = 4
          T      = Float32
          if cmp(version,"#v002") == 0 || cmp(version, "#v003")
            wdsizi = 8
            T      = Float64
          end   

#         Read the Mesh data here
          xc,yc,zc = read_re2_mesh(fid,nid0,ldimr,nelgt,wdsizi,comm)

          ncurve,curveieg,curveiside,curveparam,curvetype = read_re2_curve(fid,nid0,ldimr,nelgt,wdsizi)

          cbl,bl = read_re2_bc(fid,nid0,ldimr,nelgt,wdsizi)

          re2fld = Re2Field{T}(wdsizi,hdr,version,nelgt,ldimr,nelgv,xc,yc,zc,ncurve,curveieg,curveiside,curveparam,curvetype,cbl,bl)
         
          close(fid)
        end

        MPI.Barrier(comm)

        return re2fld 
      end     # read_re2_struct

#---------------------------------------------------------------------- 
"""
      read_re2_mesh(fid::IOStream, nid0::Int64,ldim::Int64,nelgt::Int64,wdsizi::Int64, comm::MPI.Comm)

      Read the x/y/z vertex data from the fid IOStream at MPI process ID=nid0 over the
      MPI Communicator Comm. 

      Output: xc,yc,zc

"""
      function read_re2_mesh(fid::IOStream, nid0::Int64,ldim::Int64,nelgt::Int64,wdsizi::Int64, comm::MPI.Comm)

#       Pointer to re2 data in file.
#       Header + test pattern byte length
        recpos  = 84
        seek(fid,recpos)

        rank = MPI.Comm_rank(comm)

        nc   = 2^ldim
        xc   = Array{Float64,2}(undef,nc,nelgt)
        yc   = Array{Float64,2}(undef,nc,nelgt)
        if ldim==2
          zc   = Array{Float64,2}(undef,1,1)
        else          
          zc   = Array{Float64,2}(undef,nc,nelgt)
        end  

        len     = 1 + ldim*(2^ldim)             # group + 2x4 for 2d, 3x8 for 3d
        nbytes  = len*wdsizi                    
       
        if (wdsizi == 4)
          tmp  = Vector{Float32}(undef,len)
        else
          tmp  = Vector{Float64}(undef,len)
        end
        tmp64  = Vector{Float64}(undef,len)

        for e = 1:nelgt
          if rank == nid0
            read!(fid,tmp)
            tmp64 = tmp
          end        
          if ldim == 2
            group   = tmp64[1]
#           In principle I should broadcast the data to different processors 
            xc[1:4,e] = tmp64[2:5]
            yc[1:4,e] = tmp64[6:9]
          else
#           In principle I should broadcast the data to different processors
            group   = tmp64[1]
            xc[1:8,e] = tmp64[2:9]
            yc[1:8,e] = tmp64[10:17]
            zc[1:8,e] = tmp64[18:25]
          end
        end  

        return xc,yc,zc
      end     # read_re2_mesh

#---------------------------------------------------------------------- 

      function read_re2_curve(fid::IOStream, nid0::Int64,ldim::Int64,nelt::Int64,wdsizi::Int64)

        curveieg   = Vector{Int64}(undef,1)
        curveiside = Vector{Int64}(undef,1)
        curveparam = Matrix{Float64}(undef,5,1)
        curvetype  = Vector{String}(undef,1)
        ncurve::Int64 = 0

#       Pointer to re2 data in file.
#       Header + test pattern
#        recpos  = 84
#        if wdsizi == 4
#          meshlength = (1 + ldim*(2^ldim))*4
#        else
#          meshlength = (1 + ldim*(2^ldim))*8
#        end
#        recpos = recpos + meshlength*nelt 
#        seek(fid,recpos)

        if wdsizi == 4
          nc   = read(fid,Float32)
        else
          nc   = read(fid,Float64)
        end  


        ncurve = nc
        if ncurve > 0
          curveieg   = Vector{Int64}(undef,ncurve)
          curveiside = Vector{Int64}(undef,ncurve)
          curveparam = Matrix{Float64}(undef,5,ncurve)
          curvetype  = Vector{String}(undef,ncurve)
        end 

#        len     = 2 + 1 + 5             # ieg iside curve(5) ccurve
        if (wdsizi == 4)
          tmpi  = Vector{Float32}(undef,2)
          tmpr  = Vector{Float32}(undef,5)
          tmpc  = Vector{Char}(undef,1)
        else
          tmpi  = Vector{Float64}(undef,2)
          tmpr  = Vector{Float64}(undef,5)
          tmpc  = Vector{Char}(undef,2)
        end
      
        nci = Int64(nc)
        for i in 1:nci
          
          read!(fid,tmpi)     # Read ieg, iside
          cieg          = tmpi[1]
          curveieg[i]   = Int64(cieg)
          curveiside[i] = tmpi[2]
          read!(fid,tmpr)     # Read Curve params 
          curveparam[:,i] = tmpr            
          read!(fid,tmpc)     # Read Curve type (ccurve)
          curvetype[i]    = String(tmpc)
        end  

        return ncurve,curveieg,curveiside,curveparam,curvetype 
      end     # read_re2_curve

#---------------------------------------------------------------------- 
      function read_re2_bc(fid::IOStream, nid0::Int64,ldim::Int64,nelt::Int64,wdsizi::Int64)

        nbc::Int64 = 0

        nfaces = 2*ldim
        cbl = Array{String}(undef,nfaces,nelt)
        bl  = zeros(Float64,5,nfaces,nelt)

#       Initialize Array 
        for i in 1:nelt, j in 1:nfaces
          cbl[j,i] = "E  "
        end

        nbc   = read(fid,Float64)

        if (wdsizi == 4)
          tmpi  = Vector{Int32}(undef,2)
          tmpr  = Vector{Float32}(undef,5)
          tmpc  = Vector{Char}(undef,1)
        else
          tmpi  = Vector{Float64}(undef,2)
          tmpr  = Vector{Float64}(undef,5)
          tmpc  = Vector{Char}(undef,2)
        end

#       len     = 2 + 1 + 5        # eg iside bl(5) cbl
        e::Int64  = 0 
        f::Int64  = 0
      
        bc     = Vector{Char}(undef,3)
        for i in 1:nbc
          read!(fid,tmpi)     # Read eg, iside
          e = tmpi[1]
          f = tmpi[2]
          read!(fid,tmpr)     # Read bl params 
          bl[:,f,e] = tmpr            
          read!(fid,tmpc)     # Read bl type ("E  ", etc)
          s      = String(tmpc)     # We get Length 8 for wdsizi = 8
          fill!(bc,' ')
          k      = 0
          for j in 1:length(s)
            if !isspace(s[j])
              k     = k+1
              bc[k] = s[j]
            end
          end
          cbl[f,e]  = String(bc)
        end  

#       place holder
        return cbl,bl
      end     # read_re2_bc

#---------------------------------------------------------------------- 
      function read_ma2_hdr(fid::IOStream)

        println(".ma2: Reading Header")
        
        nbytes    = 132
        hdrutf    = read(fid,nbytes)

        test         = read(fid,Float32)
        if_byte_swap = byte_swap_test(test)

        hdrC      = Char.(hdrutf)
        version   = String(hdrC[1:5])

        nilen     = 12
        i         = 5
        nelS      = String(hdrC[i+1:i+nilen])
        i         = i+nilen
        nactiveS  = String(hdrC[i+1:i+nilen])
        i         = i+nilen
        depthS    = String(hdrC[i+1:i+nilen])
        i         = i+nilen
        d2S       = String(hdrC[i+1:i+nilen])
        i         = i+nilen
        nptsS     = String(hdrC[i+1:i+nilen])
        i         = i+nilen
        nrankS    = String(hdrC[i+1:i+nilen])
        i         = i+nilen
        noutflowS = String(hdrC[i+1:i+nilen])
        i         = i+nilen

        nel       = parse(Int, nelS)
        nactive   = parse(Int, nactiveS)
        depth     = parse(Int, depthS)
        d2        = parse(Int, d2S)
        npts      = parse(Int, nptsS)
        nrank     = parse(Int, nrankS)
        noutflow  = parse(Int, noutflowS)

        hdr       = String(hdrC)

        return hdr,version,nel,nactive,depth,d2,npts,nrank,noutflow
      end     # read_ma2_hdr

#---------------------------------------------------------------------- 

      function read_ma2(f::String, nid0::Int64,comm::MPI.Comm)


        ifbcast   = false

        rank = MPI.Comm_rank(comm)

        if rank == nid0
          println("Reading $(f) on rank $rank")
        end
        

        if rank == nid0
          fid = open(f, "r")

          hdr,version,nel,nactive,depth,d2,npts,nrank,noutflow = read_ma2_hdr(fid)

          ma2hdr = ma2Hdr(version,nel,nactive,depth,d2,npts,nrank,noutflow)

          pmap, vmap = read_ma2_data(fid,nid0,nel,npts,comm)
 
          ma2data = ma2Field(pmap,vmap)
 
          close(fid)
        else
 
          ma2hdr  = ma2Hdr()
          ma2data = ma2Field()
        end  
       

        MPI.Barrier(comm)

        if ifbcast
#          buf = MPI.Buffer(hdr,132,MPI.CHAR)
#          MPI.Bcast!(buf,     nid0,comm)
          hdr       = MPI.bcast(hdr,        nid0,comm)

#          buf = MPI.Buffer(version,5,MPI.CHAR)
#          MPI.Bcast!(buf,     nid0,comm)
          version   = MPI.bcast(version,    nid0,comm)

          nel       = MPI.bcast(nel,        nid0,comm)
          nactive   = MPI.bcast(nactive,    nid0,comm)
          depth     = MPI.bcast(depth,      nid0,comm)
          d2        = MPI.bcast(d2,         nid0,comm)
          npts      = MPI.bcast(npts,       nid0,comm)
          nrank     = MPI.bcast(nrank,      nid0,comm)
          noutflow  = MPI.bcast(noutflow,   nid0,comm)
        end  

        return ma2hdr, ma2data
      end     # read_ma2

#----------------------------------------------------------------------
      function read_ma2_data(fid::IOStream, nid0::Int,nel::Int,npts::Int, comm::MPI.Comm)

#       Pointer to re2 data in file.
#       Header + test pattern byte length
        recpos  = 132+4
        seek(fid,recpos)

        rank = MPI.Comm_rank(comm)

        nvert = Int(npts/nel)
        pmap  = Vector{Int}(undef,nel)
        vmap  = Matrix{Int}(undef,nvert,nel)

        ni    = nvert+1
        line  = Vector{Int32}(undef,ni)

        for e in 1:nel
          if rank == nid0
            read!(fid,line)
            pmap[e] = line[1]
            for i in 1:nvert
              vmap[i,e] = line[i+1]
            end 
          end   # if rank == nid0
        end     # e=1:nel

        return pmap,vmap 
      end     # read_ma2_data
#----------------------------------------------------------------------

      function read_fld_struct(f::String, nid0::Int64, comm::MPI.Comm)

#       Read field file and return a NekField* Structure
#        comm = MPI.COMM_WORLD
        rank = MPI.Comm_rank(comm)

        if rank == nid0
          println("Reading $(f) on rank $rank\n")
        end
        
        MPI.Barrier(comm)

        fid = open(f, "r")

        hdr,version,wdsizi,nx,ny,nz,nel,nelgt,time,istep,fid0,nfileo,rdcode,p0th,ifprmesh,if_byte_swap = read_fld_std_hdr(fid,nid0,comm)

        buf = MPI.Buffer(hdr,length(hdr),MPI.CHAR)
        MPI.Bcast!(buf,       nid0,comm)

        buf = MPI.Buffer(version,length(version),MPI.CHAR)
        MPI.Bcast!(buf,       nid0,comm)

        wdsizi          = MPI.bcast(wdsizi,           nid0,comm)
        nx              = MPI.bcast(nx,               nid0,comm)
        ny              = MPI.bcast(ny,               nid0,comm)
        nz              = MPI.bcast(nz,               nid0,comm)
        nel             = MPI.bcast(nel,              nid0,comm)
        nelgt           = MPI.bcast(nelgt,            nid0,comm)
        time            = MPI.bcast(time,             nid0,comm)
        istep           = MPI.bcast(istep,            nid0,comm)
        fid0            = MPI.bcast(fid0,             nid0,comm)
        nfileo          = MPI.bcast(nfileo,           nid0,comm)

        buf             = MPI.Buffer(rdcode,length(rdcode),MPI.CHAR)
        MPI.Bcast!(buf,       nid0,comm)
        p0th            = MPI.bcast(p0th,             nid0,comm)
        ifprmesh        = MPI.bcast(ifprmesh,         nid0,comm)
        if_byte_swap    = MPI.bcast(if_byte_swap,     nid0,comm)

#       Read the data here        
        glnum,x,y,z,u,v,w,p,t = read_fld_data(fid, nid0,nx,ny,nz,nelgt,rdcode,wdsizi,comm)

        if wdsizi==4
#          fld = NekField8(hdr,version,wdsizi,nx,ny,nz,nel,nelgt,time,istep,fid0,nfileo,rdcode,p0th,ifprmesh,glnum,x,y,z,u,v,w,p,t)
           fld = NekField{Float32}(hdr,version,wdsizi,nx,ny,nz,nel,nelgt,time,istep,fid0,nfileo,rdcode,p0th,ifprmesh,glnum,x,y,z,u,v,w,p,t)  
         
        elseif wdsizi==8
#          fld = NekField16(hdr,version,wdsizi,nx,ny,nz,nel,nelgt,time,istep,fid0,nfileo,rdcode,p0th,ifprmesh,glnum,x,y,z,u,v,w,p,t) 
           fld = NekField{Float64}(hdr,version,wdsizi,nx,ny,nz,nel,nelgt,time,istep,fid0,nfileo,rdcode,p0th,ifprmesh,glnum,x,y,z,u,v,w,p,t)  
        else
          if (rank == nid0)
            println("Uknown word size, wdsizi=$wdsizi")
          end  
        end

        close(fid)

        return fld 
      end     # read_fld_struct

#---------------------------------------------------------------------- 

      function read_fld(f::String, nid0::Int64, comm::MPI.Comm)


        rank = MPI.Comm_rank(comm)

        if rank == nid0
          println("Reading $(f) on rank $rank)\n")
        end
        
        MPI.Barrier(comm)

        fid = open(f, "r")
        hdr,version,wdsizi,nx,ny,nz,nel,nelgt,time,istep,fid0,nfileo,rdcode,p0th,ifprmesh,if_byte_swap = read_fld_std_hdr(fid,nid0,comm)

        buf = MPI.Buffer(hdr,length(hdr),MPI.CHAR)
        MPI.Bcast!(buf,       nid0,comm)

        buf = MPI.Buffer(version,length(version),MPI.CHAR)
        MPI.Bcast!(buf,       nid0,comm)

        wdsizi          = MPI.bcast(wdsizi,           nid0,comm)
        nx              = MPI.bcast(nx,               nid0,comm)
        ny              = MPI.bcast(ny,               nid0,comm)
        nz              = MPI.bcast(nz,               nid0,comm)
        nel             = MPI.bcast(nel,              nid0,comm)
        nelgt           = MPI.bcast(nelgt,            nid0,comm)
        time            = MPI.bcast(time,             nid0,comm)
        istep           = MPI.bcast(istep,            nid0,comm)
        fid0            = MPI.bcast(fid0,             nid0,comm)
        nfileo          = MPI.bcast(nfileo,           nid0,comm)

        buf             = MPI.Buffer(rdcode,length(rdcode),MPI.CHAR)
        MPI.Bcast!(buf,       nid0,comm)
        p0th            = MPI.bcast(p0th,             nid0,comm)
        ifprmesh        = MPI.bcast(ifprmesh,         nid0,comm)
        if_byte_swap    = MPI.bcast(if_byte_swap,     nid0,comm)

#       Read the data here        
        glnum,x,y,z,u,v,w,p,t = read_fld_data(fid, nid0,nx,ny,nz,nelgt,rdcode,wdsizi,comm)

        close(fid)

        return hdr,version,wdsizi,nx,ny,nz,nel,nelgt,time,istep,fid0,nfileo,rdcode,p0th,ifprmesh,glnum,x,y,z,u,v,w,p,t
      end     # read_fld

#---------------------------------------------------------------------- 

      function read_fld_std_hdr(fid::IOStream, nid0, comm::MPI.Comm)

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

        rank = MPI.Comm_rank(comm)

        if rank == nid0
          println("Fld: Reading Header on rank $(rank)")
        
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

        test    = read(fid,Float32)

        if_byte_swap = byte_swap_test(test)

        return hdr,version,wdsize,nx,ny,nz,nel,nelgt,time,istep,fid0,nfileo,rdcode,p0th,ifprmesh,if_byte_swap
      end     # read_fld_std_hdr
#---------------------------------------------------------------------- 
      function read_fld_data(fid::IOStream, nid0::Int64,nx::Int64,ny::Int64,nz::Int64,nelgt::Int64,rdcode::String,wdsizi::Int64, comm::MPI.Comm)

#       Pointer to re2 data in file.
#       Header + test pattern byte length
        recpos  = 132+4
        seek(fid,recpos)

        rank = MPI.Comm_rank(comm)

        glnum    = Vector{Int32}(undef,nelgt)
#        glnum    = Vector{Int64}(undef,nelgt)

        read!(fid,glnum)
#        glnum = gnum

        ldim = 3
        if nz == 1
          ldim = 2
        end  

        nxyz    = nx*ny*nz
        len     = ldim*(nxyz)
        nbytes  = len*wdsizi                    
       
        if (wdsizi == 4)
          tmp   = Array{Float32}(undef,nx,ny,nz)
          tmpv  = Array{Float32}(undef,nx,ny,nz,ldim)
        else
          tmp   = Array{Float64}(undef,nx,ny,nz)
          tmpv  = Array{Float64}(undef,nx,ny,nz,ldim)
        end
        tmp64   = Array{Float64}(undef,nx,ny,nz)
        tmpv64  = Array{Float64}(undef,nx,ny,nz,ldim)

        ifxo      = false
        ifuo      = false
        ifpo      = false
        ifto      = false
        ifpso     = false
        nt        = 0
        nps       = 0

        i = 0
        for s in rdcode
          i = i+1
          if s=='X'
            ifxo = true
          end
          if s=='U'
            ifuo = true
          end
          if s=='P'
            ifpo = true
          end
          if s=='T'
            ifto = true
            nt   = 1
          end
          if s=='S'
            ifpso = true
            nps   = parse(Int64,rdcode[i+1:i+2])
          end
        end  

        nt = nt+nps
        x   = Array{Float64,4}(undef,nx,ny,nz,nelgt)
        y   = Array{Float64,4}(undef,nx,ny,nz,nelgt)
        z   = Array{Float64,4}(undef,nx,ny,nz,nelgt)

        u   = Array{Float64,4}(undef,nx,ny,nz,nelgt)
        v   = Array{Float64,4}(undef,nx,ny,nz,nelgt)
        w   = Array{Float64,4}(undef,nx,ny,nz,nelgt)

        p   = Array{Float64,4}(undef,nx,ny,nz,nelgt)
       
        t   = Array{Float64,5}(undef,nx,ny,nz,nelgt,nt)

        if (ifxo)
          for e = 1:nelgt
            if rank == nid0
              read!(fid,tmpv)
              tmpv64 = tmpv
            end
            gn = glnum[e]             
            if ldim == 2
#             In principle I should broadcast the data to different processors
              x[:,:,:,gn]  = tmpv64[:,:,:,1]
              y[:,:,:,gn]  = tmpv64[:,:,:,2]
            else
              x[:,:,:,gn]  = tmpv64[:,:,:,1]
              y[:,:,:,gn]  = tmpv64[:,:,:,2]
              z[:,:,:,gn]  = tmpv64[:,:,:,3]
            end
          end
        end  # ifxo  

        if (ifuo)
          for e = 1:nelgt
            if rank == nid0
              read!(fid,tmpv)
              tmpv64 = tmpv
            end        
            gn = glnum[e]             
            if ldim == 2
#             In principle I should broadcast the data to different processors
              u[:,:,:,gn]  = tmpv64[:,:,:,1]
              v[:,:,:,gn]  = tmpv64[:,:,:,2]
            else
              u[:,:,:,gn]  = tmpv64[:,:,:,1]
              v[:,:,:,gn]  = tmpv64[:,:,:,2]
              w[:,:,:,gn]  = tmpv64[:,:,:,3]
            end
          end
        end  # ifuo  

        if (ifpo)
          for e = 1:nelgt
            if rank == nid0
              read!(fid,tmp)
              tmp64 = tmp
            end        
#           In principle I should broadcast the data to different processors 
            gn = glnum[e]             
            p[:,:,:,gn]  = tmp64
          end
        end  # ifpo  

        if (ifto || ifpso)
          for i = 1:nt
            for e = 1:nelgt
              if rank == nid0
                read!(fid,tmp)
                tmp64 = tmp
              end        
#             In principle I should broadcast the data to different processors 
              gn = glnum[e]             
              t[:,:,:,e,i]  = tmp64
            end
          end 
        end  # ifto || ifpso 

       return glnum,x,y,z,u,v,w,p,t
      end     # read_fld_data

#---------------------------------------------------------------------- 
      function byte_swap_test(test::Float32)

        pattern::Float32  = 6.54321
        eps::Float32      = 0.00020 
        if_byte_swap      = false
         
        etest = abs(test - pattern)
        if (etest>eps) 
          if_byte_swap    = true
        end

        return if_byte_swap

      end  

#----------------------------------------------------------------------

"""
      Gen_Rema2(case::String)

      Read the case.re2 and case.ma2 files  and generate a case.rema2
      file with elements sorted according to the numbering in case.ma2.
      Boundary and Curve data also made consistent.

      Todo: Add option to output sorted.re2 and sorted.ma2

"""
      function gen_rema2(case::String, nid0::Int, comm::MPI.Comm)

        rank = MPI.Comm_rank(comm)

        re2file = case*".re2"        # String concatenation
        ma2file = case*".ma2"        # String concatenation

        hdr, map    = read_ma2(ma2file, nid0,comm)
        MPI.Barrier(comm)
        re2         = read_re2_struct(re2file,nid0,comm)

        if rank == nid0 
          nel         = re2.nelgt
          @assert hdr.nel == re2.nelgt "Total Elements don't Match. $(hdr.nel), $(nel)"

#         Get index key for sorted map          
          key         = sortperm(map.pmap)
          keyinv      = Base.copy(key)
          for i in 1:nel
            j         = key[i]
            keyinv[j] = i
          end  

          newre2      = copy(re2)
#         Renumber Coordinates and Boundary conditions        
          for i in 1:nel
            j              = key[i]
            newre2.xc[:,i] = re2.xc[:,j]
            newre2.yc[:,i] = re2.yc[:,j]
            if re2.ldimr == 3
              newre2.zc[:,i] = re2.zc[:,j]
            else
              newre2.zc[1,1] = re2.zc[1,1]
            end
            newre2.cbl[:,i]  = re2.cbl[:,j]
            newre2.bl[:,:,i] = re2.bl[:,:,j]
          end

#         Renumber Curve ieg        
          for i in 1:re2.ncurve
            el                  = re2.curveieg[i]
            newel               = keyinv[el]
            newre2.curveieg[i]  = newel
          end  

#          newmap                = ma2Field(zeros(Int64,re2.nelgt),zeros(Int,nc,nel))
          newmap                = copy(map)
          for i in 1:nel
            j                   = key[i] 
            newmap.pmap[i]      = map.pmap[j]
            newmap.vmap[:,i]    = map.vmap[:,j]
          end 

#         Generate "case.rema2.h5" file
          gen_rema2_h5(case,newre2,hdr,newmap,nid0,comm)
         

##         Write out an hdf5 file
#          fname = case*".rema2.h5"
#          fid   = h5open(fname, "w")
##         .re2 data        
#          g  = create_group(fid,"Re2")
#          g1 = create_group(g,"Params") 
#          g2 = create_group(g,"Data")        
#
#          write_dataset(g1,"ndim",newre2.ldimr)
#          write_dataset(g1,"nelgv",newre2.nelgv)
#          write_dataset(g1,"nelgt",newre2.nelgt)
#          write_dataset(g1,"wdsize",newre2.wdsize)
#
#          write_dataset(g2,"xc",newre2.xc)
#          write_dataset(g2,"yc",newre2.yc)
#          write_dataset(g2,"zc",newre2.zc)
#          write_dataset(g2,"bcs",newre2.cbl)
#          write_dataset(g2,"bcparams",newre2.bl)
#          if (newre2.ncurve>0)
#            write_dataset(g1,"ncurve",newre2.ncurve)
#            write_dataset(g2,"curveieg",newre2.curveieg)
#            write_dataset(g2,"curveiside",newre2.curveiside)
#            write_dataset(g2,"curveparam",newre2.curveparam)
#            write_dataset(g2,"curvetype",newre2.curvetype)
#          end  
#
##         .ma2 data        
#          h  = create_group(fid,"Ma2")
#          h1 = create_group(h,"Params")
#          h2 = create_group(h,"Data")
#
#          write_dataset(h1,"d2",hdr.d2)
#          write_dataset(h1,"depth",hdr.depth)
#          write_dataset(h1,"nactive",hdr.nactive)
#          write_dataset(h1,"nel",hdr.nel)
#          write_dataset(h1,"noutflow",hdr.noutflow)
#          write_dataset(h1,"npts",hdr.npts)
#          write_dataset(h1,"nrank",hdr.nrank)
#
#          write_dataset(h2,"pmap",newmap.pmap)
#          write_dataset(h2,"vmap",newmap.vmap)
#
#          close(fid)
        end       # rank == nid0  

#       For now just returning the new data        
        return nothing # newre2,newmap 
      end
#----------------------------------------------------------------------
      function gen_rema2_h5(case::String,re2::Re2Field,hdr::ma2Hdr,map::ma2Field,nid0::Int,comm::MPI.Comm)

        rank = MPI.Comm_rank(comm)

        if rank == nid0 
          nel         = re2.nelgt
          @assert hdr.nel == re2.nelgt "Total Elements don't Match. $(hdr.nel), $(nel)"

#         Write out an hdf5 file
          fname = case*".rema2.h5"
          fid   = h5open(fname, "w")
#         .re2 data        
          g  = create_group(fid,"Re2")
          g1 = create_group(g,"Params") 
          g2 = create_group(g,"Data")        

          write_dataset(g1,"ndim",re2.ldimr)
          write_dataset(g1,"nelgv",re2.nelgv)
          write_dataset(g1,"nelgt",re2.nelgt)
          write_dataset(g1,"wdsize",re2.wdsize)

          write_dataset(g2,"xc",re2.xc)
          write_dataset(g2,"yc",re2.yc)
          write_dataset(g2,"zc",re2.zc)
          write_dataset(g2,"bcs",re2.cbl)
          write_dataset(g2,"bcparams",re2.bl)
          if (re2.ncurve>0)
            write_dataset(g1,"ncurve",re2.ncurve)
            write_dataset(g2,"curveieg",re2.curveieg)
            write_dataset(g2,"curveiside",re2.curveiside)
            write_dataset(g2,"curveparam",re2.curveparam)
            write_dataset(g2,"curvetype",re2.curvetype)
          end  

#         .ma2 data        
          h  = create_group(fid,"Ma2")
          h1 = create_group(h,"Params")
          h2 = create_group(h,"Data")

          write_dataset(h1,"d2",hdr.d2)
          write_dataset(h1,"depth",hdr.depth)
          write_dataset(h1,"nactive",hdr.nactive)
          write_dataset(h1,"nel",hdr.nel)
          write_dataset(h1,"noutflow",hdr.noutflow)
          write_dataset(h1,"npts",hdr.npts)
          write_dataset(h1,"nrank",hdr.nrank)

          write_dataset(h2,"pmap",map.pmap)
          write_dataset(h2,"vmap",map.vmap)

          close(fid)
        end       # rank == nid0  

        return nothing
      end
#---------------------------------------------------------------------- 
#---------------------------------------------------------------------- 
      end   # Module JNek_IO_MPI











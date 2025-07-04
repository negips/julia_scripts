#     Author:     Prabal Negi
#     Nek based routines in SEMla

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
        ldimS   = String(hdrC[15:17])
        nelgvS  = String(hdrC[18:26])

        nelgt   = parse(Int, nelgtS)
        ldim    = parse(Int, ldimS)
        nelgv   = parse(Int, nelgvS)

        hdr     = String(hdrC)

        test    = read(fid,Float32)

        if_byte_swap = byte_swap_test(test)

        wdsize = 4
        if cmp(version,"#v002") == 0 || cmp(version, "#v003")
          wdsize = 8
        end   

        re2hdr = Re2Hdr(version,wdsize,nelgt,ldim,nelgv)    

        return hdr,re2hdr,if_byte_swap
      end     # read_re2_hdr

#----------------------------------------------------------------------
"""
      read_re2(f::String, comm::MPI.Comm, nid0::Int64)

      Read the re2 file (f) at MPI process ID=nid0 over the
      MPI Communicator Comm and output structured data Re2Field

      Output: Re2Field 

"""
      function read_re2(f::String, comm::MPI.Comm, nid0::Int64)

        re2fld = Re2Field()

        prank = MPI.Comm_rank(comm)

        if prank == nid0
          println("Reading $(f) on rank $prank")

          re2fld = read_re2(f)
          
        end

        MPI.Barrier(comm)

        return re2fld 
      end     # read_re2

#----------------------------------------------------------------------
"""
      read_re2(f::String)

      Read the re2 file (f) and output structured data Re2Field

      Output: Re2Field 

"""
      function read_re2(f::String)

        re2fld = Re2Field()

        println("Reading $(f)")

        fid = open(f, "r")
        hdr,re2hdr,if_byte_swap = read_re2_hdr(fid)

        T      = Float32
        if re2hdr.wdsize == 8
          T      = Float64
        end

#       Read the Mesh data here
        xc,yc,zc = read_re2_mesh(fid,re2hdr)

        ncurve,curveieg,curveiside,curveparam,curvetype = read_re2_curve(fid,re2hdr)

        cbl,bl = read_re2_bc(fid,re2hdr)

        re2fld = Re2Field{T}(re2hdr,xc,yc,zc,ncurve,curveieg,curveiside,curveparam,curvetype,cbl,bl)
        
        close(fid)


        return re2fld 
      end     # read_re2

#---------------------------------------------------------------------- 

"""
      read_re2_mesh(fid::IOStream, hdr::Re2Hdr)

      Read the x/y/z vertex data from the fid IOStream at MPI process ID=nid0 over the
      MPI Communicator Comm. 

      Output: xc,yc,zc

"""
      function read_re2_mesh(fid::IOStream, hdr::Re2Hdr)

#       Pointer to re2 data in file.
#       Header + test pattern byte length

        ldim      = hdr.ldim
        nelgt     = hdr.nelgt
        wdsizi    = hdr.wdsize

        recpos  = 84
        seek(fid,recpos)

        nc        = 2^ldim
        xc        = Array{Float64,2}(undef,nc,nelgt)
        yc        = Array{Float64,2}(undef,nc,nelgt)
        if ldim==2
          zc      = Array{Float64,2}(undef,1,1)
        else          
          zc      = Array{Float64,2}(undef,nc,nelgt)
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
          read!(fid,tmp)
          tmp64 = tmp
          if ldim == 2
            group   = tmp64[1]
#           In principle I should broadcast the data to different processors 
            xc[1:nc,e] = tmp64[2:5]
            yc[1:nc,e] = tmp64[6:9]
          else
#           In principle I should broadcast the data to different processors
            group      = tmp64[1]
            xc[1:nc,e] = tmp64[2:9]
            yc[1:nc,e] = tmp64[10:17]
            zc[1:nc,e] = tmp64[18:25]
          end
        end  

        return xc,yc,zc
      end     # read_re2_mesh

#---------------------------------------------------------------------- 

"""
      read_re2_curve(fid::IOStream, hdr::ReHdr)

      Read the curve data data from the fid IOStream 

      Output: ncurve,curveieg,curveiside,curveparam,curvetype 

"""
      function read_re2_curve(fid::IOStream, hdr::Re2Hdr)

        ldim       = hdr.ldim
        nelt       = hdr.nelgt
        wdsizi     = hdr.wdsize

        curveieg   = Vector{Int64}(undef,1)
        curveiside = Vector{Int64}(undef,1)
        curveparam = Matrix{Float64}(undef,5,1)
        curvetype  = Vector{String}(undef,1)
        ncurve     = 0

        if wdsizi == 4
          T    = Float32
        else
          T    = Float64
        end  

        nc   = read(fid,T)

        ncurve = Int(nc)
        if ncurve > 0
          curveieg   = Vector{Int64}(undef,ncurve)
          curveiside = Vector{Int64}(undef,ncurve)
          curveparam = Matrix{Float64}(undef,5,ncurve)
          curvetype  = Vector{String}(undef,ncurve)
        end 

#        len     = 2 + 1 + 5             # ieg iside curve(5) ccurve
        tmpi  = Vector{T}(undef,2)
        tmpr  = Vector{T}(undef,5)
       
        if (T == Float32)
          tmpc  = Vector{Char}(undef,1)
        else
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
      function read_re2_bc(fid::IOStream, hdr::Re2Hdr)

        ldim      = hdr.ldim
        nelt      = hdr.nelgt
        wdsizi    = hdr.wdsize

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
"""
      read_ma2(f::String, comm::MPI.Comm, nid0::Int64)

      Read the ma2 file (f) at MPI process ID=nid0 over the
      MPI Communicator Comm.

      Output: Ma2fld

"""
      function read_ma2(f::String, comm::MPI.Comm, nid0::Int64)

        prank = MPI.Comm_rank(comm)

        ma2fld = Ma2Field()

        if prank == nid0
          println("Reading $(f) on rank $prank")

          ma2fld  = read_ma2(f)
 
        end  

        MPI.Barrier(comm)

        return ma2fld
      end     # read_ma2

#----------------------------------------------------------------------
"""
      read_ma2(f::String)

      Read the ma2 file (f) and output the Ma2Field

      Output: Ma2fld

"""
      function read_ma2(f::String)

        println("Reading $(f)")

        fid = open(f, "r")

        hdr,ma2hdr = read_ma2_hdr(fid)

        pmap, vmap = read_ma2_data(fid,ma2hdr)
 
        ma2fld     = Ma2Field(ma2hdr,pmap,vmap)
 
        close(fid)

        return ma2fld
      end     # read_ma2

#----------------------------------------------------------------------

"""
      read_ma2_hdr(fid::IOStream)

      Read the header of IOStream from a .ma2 file.

"""
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

        ma2hdr    = Ma2Hdr(version,nel,nactive,depth,d2,npts,nrank,noutflow)

        return hdr,ma2hdr
      end     # read_ma2_hdr

#---------------------------------------------------------------------- 

"""
      read_ma2_data(fid::IOStream, hdr::Ma2Hdr)

      Read the Element and vertex map from the fid IOStream 

      Output: pmap,vmap

"""
      function read_ma2_data(fid::IOStream, hdr::Ma2Hdr)

        nel       = hdr.nel
        npts      = hdr.npts

#       Pointer to re2 data in file.
#       Header + test pattern byte length
        recpos  = 132+4
        seek(fid,recpos)

        nvert = Int(npts/nel)
        pmap  = Vector{Int}(undef,nel)
        vmap  = Matrix{Int}(undef,nvert,nel)

        ni    = nvert+1
        line  = Vector{Int32}(undef,ni)

        for e in 1:nel
          read!(fid,line)
          pmap[e] = line[1]
          for i in 1:nvert
            vmap[i,e] = line[i+1]
          end 
        end     # e=1:nel

        return pmap,vmap 
      end     # read_ma2_data
#----------------------------------------------------------------------
"""
      read_fld(f::String, comm::MPI.Comm, nid0::Int64)

      Read the Nek Field file (f) at MPI process ID=nid0 over the
      MPI Communicator Comm and output structured data NekField

      Output: NekField 

"""
      function read_fld(f::String, comm::MPI.Comm, nid0::Int64)

#       Read field file and return a NekField Structure
        prank = MPI.Comm_rank(comm)

        T   = Float32
        fld = NekField(T)

        MPI.Barrier(comm)

        if prank == nid0
          println("Reading $(f) on rank $prank\n")
      
          fld = read_fld(f)

        end

        return fld 
      end     # read_fld

#---------------------------------------------------------------------- 
"""
      read_fld(f::String)

      Read the Nek Field file (f) and output structured data NekField

      Output: NekField 

"""
      function read_fld(f::String)

#       Read field file and return a NekField Structure

        println("Reading $(f)\n")

        fid = open(f, "r")

        # Read Header
        hdrS,fldhdr,if_byte_swap = read_fld_std_hdr(fid)
        T = Float32
        if fldhdr.wdsize == 8
          T = Float64
        end  
         
#       Read data
        glnum,x,y,z,u,v,w,p,t = read_fld_data(fid,fldhdr)
        fld = NekField{T}(fldhdr,glnum,x,y,z,u,v,w,p,t)  

        return fld 
      end     # read_fld

#---------------------------------------------------------------------- 

      function read_fld_std_hdr(fid::IOStream)
        

        println("Fld: Reading Header")
        
        nbytes        = 132
        hdrutf        = read(fid,nbytes)
        test          = read(fid,Float32)
        if_byte_swap  = byte_swap_test(test)

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

        T      = Float32
        if wdsize==8
          T    = Float64
        end  

        fldhdr = NekFldHdr{T}(version,wdsize,nx,ny,nz,nel,nelgt,time,istep,fid0,nfileo,rdcode,p0th,ifprmesh)

        return hdr,fldhdr,if_byte_swap
      end     # read_fld_std_hdr
#---------------------------------------------------------------------- 
      function read_fld_data(fid::IOStream, fldhdr::NekFldHdr)

        nx      = fldhdr.lx1
        ny      = fldhdr.ly1
        nz      = fldhdr.lz1
        nelgt   = fldhdr.nelgt
        rdcode  = fldhdr.rdcode
        wdsizi  = fldhdr.wdsize

#       Pointer to re2 data in file.
#       Header + test pattern byte length
        recpos  = 132+4
        seek(fid,recpos)

        glnum    = Vector{Int32}(undef,nelgt)

        read!(fid,glnum)

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
            read!(fid,tmpv)
            tmpv64 = tmpv
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
            read!(fid,tmpv)
            tmpv64 = tmpv
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
            read!(fid,tmp)
            tmp64 = tmp
#           In principle I should broadcast the data to different processors 
            gn = glnum[e]             
            p[:,:,:,gn]  = tmp64
          end
        end  # ifpo  

        if (ifto || ifpso)
          for s = 1:nt
            for e = 1:nelgt
              read!(fid,tmp)
              tmp64 = tmp
#             In principle I should broadcast the data to different processors 
              gn = glnum[e]             
              t[:,:,:,e,s]  = tmp64
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


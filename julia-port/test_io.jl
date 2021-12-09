# testing the re2/map/par readers
#

      using Revise

#      include("JNek_PARALLEL.jl")         # JNek_PARALLEL
      include("JNek_IO.jl")               # JNek_IO

#      import JNek_IO
      
      using MPI

      nid0  = 0

      if (MPI.Initialized() == false)
        MPI.Init()
      end  
        
      comm = MPI.COMM_WORLD
      rank = MPI.Comm_rank(comm)

      re2   = "channelp.re2"
      map   = "channelp.ma2"
      fld   = "taylor0.f00001"

      hdr,version,nelgt,ldimr,nelgv,xc,yc,ncurve,curveieg,curveiside,curveparam,curvetype,cbl,bl = JNek_IO.read_re2(re2,nid0)

      hdr,version,wdsize,nx,ny,nz,nel,nelgt,time,istep,fid0,nfileo,rdcode,p0th,ifprmesh,glnum,x,y,z,u,v,w,p,t = JNek_IO.read_fld(fld,MPI,nid0)
      if rank == 0
        println("Done")
      end  

#      MPI.Finalize()

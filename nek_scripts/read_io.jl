# testing the re2/map/par readers
#

      using Revise
      using PyCall
      using PyPlot 

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

      rea   = "box.rea"
      re2   = "test.re2"
      map   = "channelp.ma2"
#      fld   = "inttaylor0.f00001"
      fld1   = "tstlevelset0.f00001"
      fld2   = "newlevelset0.f00001"

      wdsizi,hdr,version,nelgt,ldimr,nelgv,xc,yc,ncurve,curveieg,curveiside,curveparam,curvetype,cbl,bl = JNek_IO.read_re2(re2,nid0)

#      hdr,version,wdsize,nx,ny,nz,nel,nelgt,time,istep,fid0,nfileo,rdcode,p0th,ifprmesh,glnum,x1,y1,z1,u1,v1,w1,p1,t1 = JNek_IO.read_fld(fld1,MPI,nid0)

#      hdr,version,wdsize,nx,ny,nz,nel,nelgt,time,istep,fid0,nfileo,rdcode,p0th,ifprmesh,glnum,x2,y2,z2,u2,v2,w2,p2,t2 = JNek_IO.read_fld(fld2,MPI,nid0)

#      du = u1[:,:,:,1] .- u2[:,:,:,1]
#      uv = du[1,:,1,1]
#      y1 = y1[1,:,1,1];
   
      x1 = xc[:]
      y1 = yc[:]
      scatter(x1,y1)
      if rank == 0
        println("Done")
      end  

#      MPI.Finalize()

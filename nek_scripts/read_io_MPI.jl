# testing the re2/map/par readers
#

      using Revise
      using PyCall
      using PyPlot 
      
      include("JNek_IO_MPI.jl")            # JNek_IO
      
      using MPI

      nid0  = 0

      if !MPI.Initialized()
        MPI.Init()
      end  
        
      comm    = MPI.COMM_WORLD
      rank    = MPI.Comm_rank(comm)

      rea     = "box.rea"
      re2     = "three.re2"
      ma2     = "three.ma2"
      file1   = "cyl0.f00001"

#      wdsizi,hdr,version,nelgt,ldimr,nelgv,xc,yc,zc,ncurve,curveieg,curveiside,curveparam,curvetype,cbl,bl = JNek_IO.read_re2(re2,nid0)
      re2fld1 = JNek_IO.read_re2_struct(re2,nid0,comm)

      ma2hdr,ma2data  = JNek_IO.read_ma2(ma2,nid0,comm)

      fld  = JNek_IO.read_fld_struct(file1,MPI,nid0, comm)

#      hdr,version,wdsize,nx,ny,nz,nel,nelgt,time,istep,fid0,nfileo,rdcode,p0th,ifprmesh,glnum,x2,y2,z2,u2,v2,w2,p2,t2 = JNek_IO.read_fld(file12,MPI,nid0)
#       e    = 1; 
#       varx = fld.x[:,:,1,e]; 
#       vary = fld.y[:,:,1,e]; 
#       varu = fld.u[:,:,1,e]; 
#       varv = fld.v[:,:,1,e];
#       varp = fld.p[:,:,1,e];
#       vart = fld.t[:,:,1,e,1];
#
#       for i in [1 4]
#         xe   = varx[:,i];
#         pe   = varp[:,i];
#         ve   = varv[:,i];
##         plot(xe,ve)
#       end  
#
#      if rank == 0
#        println("Done")
#      end  

#      MPI.Finalize()

#      MPI.Initialized()

      println("Done. MPI not Finalized")


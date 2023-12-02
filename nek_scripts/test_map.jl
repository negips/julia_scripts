# testing the re2/map/par readers
#

      using Revise
      using PyCall
      using PyPlot 
      
      include("JNek_IO.jl")            # JNek_IO
      include("NekTools.jl")
      include("$(JULIACOMMON)/MoveFigure.jl")
      
      using MPI

      close("all")

      nid0  = 0

      if !MPI.Initialized()
        MPI.Init()
      end  
        
      comm    = MPI.COMM_WORLD
      rank    = MPI.Comm_rank(comm)

      rea     = "box.rea"
      re2     = "cyl.re2"
      ma2     = "cyl.ma2"
      file1   = "cyl0.f00001"

      re2fld          = JNek_IO.read_re2_struct(re2,nid0,comm)
      ma2hdr,ma2data  = JNek_IO.read_ma2(ma2,nid0,comm)

      xc = zeros(Float64,size(re2fld.xc))
      yc = similar(xc)

      nel   = re2fld.nelgv
      ldim  = re2fld.ldimr
      nc    = 2^ldim

      fld  = JNek_IO.read_fld_struct(file1,nid0, comm)
      
      for e in 1:nel
        for i in 1:nc
          j       = NekPrepToSymm(i)
          xc[j,e] = re2fld.xc[i,e]
          yc[j,e] = re2fld.yc[i,e]
        end
      end


      xc2 = copy(re2fld.xc)
      yc2 = copy(re2fld.yc)

      cm2   = get_cmap("RdBu_r");   #_r for reversed.
      h1    = figure(num=1)

#      pcm   = pcolormesh(xc,yc,Pm)
#      pcm   = surf(Xc,Yc,Pm)
#      pcm.set_cmap(cm2)
      u0    = minimum(fld.u[:])
      u1    = maximum(fld.u[:])

#      matplotlib.colors.Normalize(vmin=u0, vmax=u1)

      lx1 = fld.lx1
      ly1 = fld.ly1

      x2d  = zeros(Float64,2,2)
      y2d  = zeros(Float64,2,2)
      f2d  = zeros(Int64,2,2)

      c1   = CartesianIndex(1,1)
      c2   = CartesianIndex(lx1,1)
      c3   = CartesianIndex(1,ly1)
      c4   = CartesianIndex(lx1,ly1)

      p0    = minimum(ma2data.pmap)
      p1    = maximum(ma2data.pmap)


      h1 = figure(num=1)
      MoveFigure(h1,1750,10)

      for e in 1:nel
        x2d[1,1] = fld.x[c1,1,e]
        x2d[2,1] = fld.x[c2,1,e]
        x2d[1,2] = fld.x[c3,1,e]
        x2d[2,2] = fld.x[c4,1,e]

        y2d[1,1] = fld.y[c1,1,e]
        y2d[2,1] = fld.y[c2,1,e]
        y2d[1,2] = fld.y[c3,1,e]
        y2d[2,2] = fld.y[c4,1,e]

        pm  = ma2data.pmap[e]
        fill!(f2d,pm)

        pcm   = pcolormesh(x2d,y2d,f2d,vmin=p0,vmax=p1)

#        surf(x2d,y2d,f2d,vmin=u0,vmax=u1)
#        pcm   = pcolormesh(x2d,y2d,f2d)
        pcm.set_cmap(cm2)
      end  


      h2 = figure(num=2)
      MoveFigure(h2,1750,800)

      newre2,newmap = JNek_IO.gen_rema2("cyl",nid0,comm)

      for e in 1:nel
        x2d[1,1] = newre2.xc[1,e]
        x2d[2,1] = newre2.xc[2,e]
        x2d[1,2] = newre2.xc[4,e]
        x2d[2,2] = newre2.xc[3,e]

        y2d[1,1] = newre2.yc[1,e]
        y2d[2,1] = newre2.yc[2,e]
        y2d[1,2] = newre2.yc[4,e]
        y2d[2,2] = newre2.yc[3,e]

        pm2  = newmap.pmap[e]
        fill!(f2d,pm2)

        pcm2   = pcolormesh(x2d,y2d,f2d,vmin=p0,vmax=p1)

#        surf(x2d,y2d,f2d,vmin=u0,vmax=u1)
#        pcm   = pcolormesh(x2d,y2d,f2d)
        pcm2.set_cmap(cm2)
      end  





#      if rank == 0
#        println("Done")
#      end  

#      MPI.Finalize()

#      MPI.Initialized()

      if rank == 0
        println("Done. MPI not Finalized")
      end  


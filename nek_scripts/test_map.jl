# testing the re2/map/par readers
#

      using Revise
      using PyCall
      using PyPlot 
      
      include("JNek_IO_MPI.jl")            # JNek_IO
      include("NekTools.jl")
      
      using MPI

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

      for e in 1:nel
        for i in 1:nc
          j       = NekPrepToSymm(i)
          xc[j,e] = re2fld.xc[i,e]
          yc[j,e] = re2fld.yc[i,e]
        end
      end


#      fld  = JNek_IO.read_fld_struct(file1,MPI,nid0, comm)

#      if rank == 0
#        println("Done")
#      end  

#      MPI.Finalize()

#      MPI.Initialized()

      if rank == 0
        println("Done. MPI not Finalized")
      end  


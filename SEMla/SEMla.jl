#     Author:     Prabal Negi
#

      module SEMla

      using MPI
      using HDF5

#     Include for basic definitions: abstract types, structures, constructors, extensions      
      include("SEMla_Abstract.jl")        # 
      include("SEMla_Structs.jl")         # 
      include("SEMla_Constructors.jl")    #
      include("SEMla_Extends.jl")         # native function extensions

#     Include files for defined functions
      include("SEMla_Nek.jl")             # Exports the Nek related functions


#     Abstraact Definitions      
      export AbstractTensorfield,
             AbstractIsoParTensorField


      export NekField,
             Re2Field,
             TwoTensorField


      export gen_rema2,
             redistribute_mesh

#     Nek based functions      
      export read_re2_hdr,
             read_re2,
             read_re2_struct,
             read_ma2_hdr,
             read_ma2,
             read_fld,
             read_fld_struct

#----------------------------------------------------------------------  

      function __init__()

        return nothing
      end 

#----------------------------------------------------------------------

      function init()

#       Can't reinitialize MPI        

        return nothing
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
          keyinv      = Base.copy(key)          # Inverse key
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

          newmap                = copy(map)
          for i in 1:nel
            j                   = key[i] 
            newmap.pmap[i]      = map.pmap[j]
            newmap.vmap[:,i]    = map.vmap[:,j]
          end 

#         Generate "case.rema2.h5" file
          gen_rema2_h5(case,newre2,hdr,newmap,nid0,comm)

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
      end         # gen_rema2_h5
#---------------------------------------------------------------------- 
"""
      Distribute_Mesh(case::String)

      Read the case.re2 and case.ma2 files  and generate a case.rema2
      file with elements sorted according to the numbering in case.ma2.
      Boundary and Curve data also made consistent.

"""
      function Distribute_Mesh(casename::String,nid0::Int,comm::MPI.Comm)

        rank        = MPI.Comm_rank(comm)
        np          = MPI.Comm_size(comm)

        h5name      = casename*".rema2.h5"
        gnel        = 0
        ndim        = 0
        wdsize      = 4
        T           = Float32

        if rank == nid0
          fid      = h5open(h5name, "r")

          g        = fid["Re2"]
          g1       = g["Params"]
          g2       = g["Data"]

          h        = fid["Ma2"]
          h1       = h["Params"]
          h2       = h["Data"]

          gnel     = read(g1,"nelgt")
          ndim     = read(g1,"ndim")
          wdsize   = read(g1,"wdsize")

        end  

        gnel       = MPI.bcast(gnel, nid0, comm)
        ndim       = MPI.bcast(ndim, nid0, comm)
        wdsize     = MPI.bcast(wdsize, nid0, comm)

        if wdsize == 8
          T = Float64
        end  

        # Remaining Elements
        rem   = mod(gnel,np)  
        nel2  = gnel - rem
        lnel0 = floor(Int,nel2/np)
        # Add remaining elements to the end ranks
        last_rank = np-1     
        pextras = np - rem 
   
        lnel  = lnel0
        if rank >= pextras
          lnel = lnel0 + 1
        end

        # Allocate memory for xc,yc
        nc    = 2^ndim
        xc    = zeros(T,nc,lnel)
        yc    = zeros(T,nc,lnel)

        # Temporary array to store data for other processors        
        xcb   = zeros(T,nc,lnel0+1)
        ycb   = zeros(T,nc,lnel0+1)

        if rank == nid0
          # Read from the case.rema2.h5 file        
          xcg  = read(g2,"xc")
          ycg  = read(g2,"yc")
         
          j    = 0
          nid  = 0
          jend = 0
          while nid < np
            # global jend, j, nid
            jend  = jend + lnel0
            if nid >= pextras
              jend   = jend + 1
            end

            i    = 0
            while j < jend
              j = j + 1
              i = i + 1

              ind1  = CartesianIndices((1:nc,i:i))
              ind2  = CartesianIndices((1:nc,j:j)) 
              copyto!(xcb,ind1,xcg,ind2)
              copyto!(ycb,ind1,ycg,ind2)

            end     # while j < jend

            if nid == nid0
              ind1  = CartesianIndices((1:nc,1:lnel))
              copyto!(xc,ind1,xcb,ind1)
              copyto!(yc,ind1,ycb,ind1)
            else  
              destn        = nid
              send_status1 = MPI.Isend(xcb,comm;dest=destn)

              destn        = nid
              send_status2 = MPI.Isend(ycb,comm;dest=destn)
            end     # if nid == nid0
            nid = nid + 1
          end         # while nid < np 

        else  # rank != nid0
          recv_status1 = MPI.Irecv!(xcb,comm; source=nid0)
          recv_status2 = MPI.Irecv!(ycb,comm; source=nid0)

          stat1 =  MPI.Wait!(recv_status1)
          stat2 =  MPI.Wait!(recv_status2)

          ind3  = CartesianIndices((1:nc,1:lnel))
          copyto!(xc,ind3,xcb,ind3)
          copyto!(yc,ind3,ycb,ind3)
        end   # if rank == nid0 

        MPI.Barrier(comm)

#       Checks!        
#       This should be commented out after checks        
#       Write out an hdf5 file
        ftest = casename*"_p$rank.h5"
        fid2  = h5open(ftest, "w")
        gg    = create_group(fid2,"Coords")
        write_dataset(gg,"xc",xc)
        write_dataset(gg,"yc",yc)
        close(fid2)

        println("\n Writing Nel = $lnel, Elements, on Rank=$rank to $ftest\n")

        return nothing
      end  # Distribute_Mesh 
#----------------------------------------------------------------------  


#---------------------------------------------------------------------- 
      end   # Module SEMla











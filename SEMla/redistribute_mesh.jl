# testing the re2/map/par readers

      function redistribute_mesh(casename::String,nid0::Int,comm::MPI.Comm)


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
#            global jend, j, nid
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

#       Write out an hdf5 file
        ftest = casename*"_p$rank.h5"
        fid2  = h5open(ftest, "w")
        gg    = create_group(fid2,"Coords")
        write_dataset(gg,"xc",xc)
        write_dataset(gg,"yc",yc)
        close(fid2)

        println("\n Writing Nel = $lnel, Elements, on Rank=$rank to $ftest\n")

      end        

#-------------------------------------------------- 



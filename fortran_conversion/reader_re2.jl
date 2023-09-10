using FortranFiles
using OffsetArrays
using Parameters
using Printf

# COMMON /CTMP1/ VI
@static if !@isdefined(ctmp1_COMMON_STRUCT)
@with_kw mutable struct ctmp1_COMMON_STRUCT
    vi::Matrix{Int32}
end
end

#ctmp1 = ctmp1_COMMON_STRUCT()


# COMMON /SCRNS/ BUFR
@static if !@isdefined(scrns_COMMON_STRUCT)
@with_kw mutable struct scrns_COMMON_STRUCT
    bufr::Matrix{Int32}
end
end

#scrns = scrns_COMMON_STRUCT()


# COMMON /NEKMPI/ NIDD,NPP,NEKCOMM,NEKGROUP,NEKREAL
@static if !@isdefined(nekmpi_COMMON_STRUCT)
@with_kw mutable struct nekmpi_COMMON_STRUCT
    nidd::Int64
    npp::Int64
    nekcomm::Int64
    nekgroup::Int64
    nekreal::Int64
end
end

#nekmpi = nekmpi_COMMON_STRUCT()


#-----------------------------------------------------------------------
      function read_redata(2)(ifbswap)  # .re2 reader

      include("SIZE")
      include("TOTAL")
      include("RESTART")
      include("CTIMER")

#      logical ifbswap
#      integer idummy[100]

#      common /nekmpi/ nidd,npp,nekcomm,nekgroup,nekreal
      global nekmpi
      @unpack NIDD, NPP, NEKCOMM, NEKGROUP, NEKREAL = nekmpi


      etime0 = dnekclock_sync()

                  ibc = 2
      if (ifflow) ibc = 1 end

                  nfldt = 1
      if (ifheat) nfldt = 2+npscal end
      if (ifmhd) nfldt = 2+npscal+1 end

      # first field to read
      if (param(33) > 0) ibc = trunc(Int, param(33)) end

      # number of fields to read
      if (param(32) > 0) nfldt = ibc + trunc(Int, param(32)) - 1 end

      blank(cbc, 3*SIZE(cbc))
      rzero(bc , SIZE(bc))

#ifndef NOMPIIO
      fgslib_crystal_setup(cr_re2, nekcomm, np)

      byte_open_mpi(re2fle, fh_re2, true, ierr)
      err_chk(ierr, " Cannot open .re2 file!\$")

      readp_remesh(2)(ifbswap)
      readp_recurve(2)(ifbswap)
      for ifield = ibc:nfldt
         readp_rebc(2)(cbc(1, 1, ifield), bc(1, 1, 1, ifield), ifbswap)
      end

      fgslib_crystal_free(cr_re2)
      byte_close_mpi(fh_re2, ierr)
#else
      byte_open(re2fle, ierr)
      byte_read(idummy, 21, ierr) # skip hdr+endian code 

      bin_rdmesh(1)(ifbswap)
      bin_rdcurve(1)(ifbswap)
      for ifield = ibc:nfldt
         bin_rdbc(1)(cbc(1, 1, ifield), bc(1, 1, 1, ifield), ifbswap)
      end

      byte_close(ierr)
#endif

      etime_t = dnekclock_sync() - etime0
      if (nio == 0) @printf(stdout, "%sg9.2%s\n\n", " done :: read .re2 file   ",
                         etime_t, " sec")
      end


      @label Lreturn
      @pack! NEKMPI = NIDD, NPP, NEKCOMM, NEKGROUP, NEKREAL
      return nothing
      end
#-----------------------------------------------------------------------
      function readp_remesh(2)(ifbswap) # version 2 of .re2 reader

      include("SIZE")
      include("TOTAL")

#      logical ifbswap

                nrmax = lelt             # maximum number of records
                lrs   = 1+ldim*Complex(2^ldim) # record size: group x(:,c) ...
                li    = 2*lrs+2

#      integer         bufr[li-2,nrmax]
#      common /scrns/  bufr
      global scrns
      @unpack BUFR = scrns

#      integer         vi[li  ,nrmax]
#      common /ctmp1/  vi
      global ctmp1
      @unpack VI = ctmp1

#      integer*8       lre2off_b,dtmp8
#      integer*8       nrg

      if (nio == 0) println(stdout, " preading mesh ") end

      nrg       = nelgt
      nr        = nelt
      irankoff  = igl_running_sum(nr) - nr
      dtmp8     = irankoff
      re2off_b  = 84 # set initial offset (hdr + endian)
      lre2off_b = re2off_b + dtmp8*lrs*wdsizi
      lrs4      = lrs*wdsizi/4

      # read coordinates from file
      nwds4r = nr*lrs4
      byte_set_view(lre2off_b, fh_re2)
      byte_read_mpi(bufr, nwds4r, -1, fh_re2, ierr)
      re2off_b = re2off_b + nrg*4*lrs4
      if (ierr > 0) @goto L100 end

      # pack buffer
      for i = 1:nr
         jj      = (i-1)*lrs4 + 1
         ielg    = irankoff + i # elements are stored in global order
         vi[1,i] = gllnid(ielg)
         vi[2,i] = ielg
         icopy(vi[3,i], bufr[jj,1], lrs4)
      end

      # crystal route nr real items of size lrs to rank vi(key,1:nr)
      n   = nr
      key = 1
      fgslib_crystal_tuple_transfer(cr_re2, n, nrmax, vi, li, vl, 0, vr, 0,
                                         key)

      # unpack buffer
      ierr = 0
      if n > nrmax
         ierr = 1
          @goto L100
      end

      for i = 1:n
         iel = gllel(vi[2,i])
         icopy(bufr, vi[3,i], lrs4)
         buf_to_xyz(bufr, iel, ifbswap, ierr)
      end

      @label L100
      err_chk(ierr, "Error reading .re2 mesh\$")


      @label Lreturn
      return nothing
      end
#-----------------------------------------------------------------------
      function readp_recurve(2)(ifbswap)

      include("SIZE")
      include("TOTAL")

#      logical ifbswap

#      common /nekmpi/ nidd,npp,nekcomm,nekgroup,nekreal
      global nekmpi
      @unpack NIDD, NPP, NEKCOMM, NEKGROUP, NEKREAL = nekmpi

                nrmax = 12*lelt # maximum number of records
                lrs   = 2+1+5   # record size: eg iside curve(5) ccurve
                li    = 2*lrs+1

#      integer         bufr[li-1,nrmax]
#      common /scrns/  bufr
      global scrns
      @unpack BUFR = scrns

#      integer         vi[li  ,nrmax]
#      common /ctmp1/  vi
      global ctmp1
      @unpack VI = ctmp1

#      integer*8       lre2off_b,dtmp8
#      integer*8       nrg
#      integer*4       nrg4[2]


      # read total number of records
      nwds4r    = 1*wdsizi/4
      lre2off_b = re2off_b
      byte_set_view(lre2off_b, fh_re2)
      byte_read_mpi(nrg4, nwds4r, -1, fh_re2, ierr)
      if (ierr > 0) @goto L100 end

      if wdsizi == 8
         if (ifbswap) byte_reverse8(nrg4, nwds4r, ierr) end
         COPY(dnrg, nrg4, 1)
         nrg = dnrg
      else
         if (ifbswap) byte_reverse(nrg4, nwds4r, ierr) end
         nrg = nrg4[1]
      end
      re2off_b = re2off_b + 4*nwds4r

      if (nrg == 0) @goto Lreturn end
      if (nio == 0) println(stdout, " preading curved sides ") end

      # read data from file
      dtmp8 = np
      nr = nrg/dtmp8
      for i = 0:mod(nrg, dtmp8)-1
         if (i == nid) nr = nr + 1 end
      end
      irankoff  = igl_running_sum(nr) - nr
      dtmp8     = irankoff
      lre2off_b = re2off_b + dtmp8*lrs*wdsizi
      lrs4      = lrs*wdsizi/4

      nwds4r = nr*lrs4
      byte_set_view(lre2off_b, fh_re2)
      byte_read_mpi(bufr, nwds4r, -1, fh_re2, ierr)

      re2off_b = re2off_b + nrg*4*lrs4
      if (ierr > 0) @goto L100 end

      # pack buffer
      for i = 1:nr
         jj = (i-1)*lrs4 + 1

         if ifbswap
           lrs4s = lrs4 - wdsizi/4 # words to swap (last is char)
           if (wdsizi == 8) byte_reverse8(bufr[jj,1], lrs4s, ierr) end
           if (wdsizi == 4) byte_reverse(bufr[jj,1], lrs4s, ierr) end
         end

         ielg = bufr[jj,1]
         if (wdsizi == 8) copyi4(ielg, bufr[jj,1], 1) end

         if (ielg <= 0 || ielg > nelgt) @goto L100 end
         vi[1,i] = gllnid(ielg)

         icopy(vi[2,i], bufr[jj,1], lrs4)
      end

      # crystal route nr real items of size lrs to rank vi(key,1:nr)
      n    = nr
      key  = 1
      fgslib_crystal_tuple_transfer(cr_re2, n, nrmax, vi, li, vl, 0, vr, 0,
                                         key)

      # unpack buffer
      if (n > nrmax) @goto L100 end
      for i = 1:n
         icopy(bufr, vi[2,i], lrs4)
         buf_to_curve(bufr)
      end

      @goto Lreturn

      @label L100
      ierr = 1
      err_chk(ierr, "Error reading .re2 curved data\$")


      @label Lreturn
      @pack! NEKMPI = NIDD, NPP, NEKCOMM, NEKGROUP, NEKREAL
      return nothing
      end
#-----------------------------------------------------------------------
      function readp_rebc(2)(cbl, bl, ifbswap)

      include("SIZE")
      include("TOTAL")

#      character*3  cbl[  6,lelt]
#      REAL         bl[5,6,lelt]
#      logical      ifbswap

                nrmax = 6*lelt # maximum number of records
                lrs   = 2+1+5  # record size: eg iside bl(5) cbl
                li    = 2*lrs+1

#      integer         bufr[li-1,nrmax]
#      common /scrns/  bufr
      global scrns
      @unpack BUFR = scrns

#      integer         vi[li  ,nrmax]
#      common /ctmp1/  vi
      global ctmp1
      @unpack VI = ctmp1

#      integer*8       lre2off_b,dtmp8
#      integer*8       nrg
#      integer*4       nrg4[2]


      # fill up with default
      for iel = 1:nelt
      for k = 1:6
         cbl[k,iel] = "E  "
      end
      end

      # read total number of records
      nwds4r    = 1*wdsizi/4
      lre2off_b = re2off_b
      byte_set_view(lre2off_b, fh_re2)
      byte_read_mpi(nrg4, nwds4r, -1, fh_re2, ierr)
      if (ierr > 0) @goto L100 end

      if wdsizi == 8
         if (ifbswap) byte_reverse8(nrg4, nwds4r, ierr) end
         COPY(dnrg, nrg4, 1)
         nrg = dnrg
      else
         if (ifbswap) byte_reverse(nrg4, nwds4r, ierr) end
         nrg = nrg4[1]
      end
      re2off_b = re2off_b + 4*nwds4r

      if (nrg == 0) @goto Lreturn end
      if (nio == 0) println(stdout, " preading bc for ifld", ifield) end

      # read data from file
      dtmp8 = np
      nr = nrg/dtmp8
      for i = 0:mod(nrg, dtmp8)-1
         if (i == nid) nr = nr + 1 end
      end
      irankoff  = igl_running_sum(nr) - nr
      dtmp8     = irankoff
      lre2off_b = re2off_b + dtmp8*lrs*wdsizi
      lrs4      = lrs*wdsizi/4

      nwds4r = nr*lrs4
      byte_set_view(lre2off_b, fh_re2)
      byte_read_mpi(bufr, nwds4r, -1, fh_re2, ierr)

      re2off_b = re2off_b + nrg*4*lrs4
      if (ierr > 0) @goto L100 end

      # pack buffer
      for i = 1:nr
         jj = (i-1)*lrs4 + 1

         if ifbswap
           lrs4s = lrs4 - wdsizi/4 # words to swap (last is char)
           if (wdsizi == 8) byte_reverse8(bufr[jj,1], lrs4s, ierr) end
           if (wdsizi == 4) byte_reverse(bufr[jj,1], lrs4s, ierr) end
         end

         ielg = bufr[jj,1]
         if (wdsizi == 8) copyi4(ielg, bufr[jj,1], 1) end

         if (ielg <= 0 || ielg > nelgt) @goto L100 end
         vi[1,i] = gllnid(ielg)

         icopy(vi[2,i], bufr[jj,1], lrs4)
      end

      # crystal route nr real items of size lrs to rank vi(key,1:nr)
      n    = nr
      key  = 1
      fgslib_crystal_tuple_transfer(cr_re2, n, nrmax, vi, li, vl, 0, vr, 0,
                                         key)

      # unpack buffer
      if (n > nrmax) @goto L100 end
      for i = 1:n
         icopy(bufr, vi[2,i], lrs4)
         buf_to_bc(cbl, bl, bufr)
      end

      @goto Lreturn

      @label L100
      ierr = 1
      err_chk(ierr, "Error reading .re2 boundary data\$")


      @label Lreturn
      return nothing
      end
#-----------------------------------------------------------------------
      function buf_to_xyz(buf, e, ifbswap, ierr)# version 1 of binary reader

      include("SIZE")
      include("TOTAL")
#      logical ifbswap

#      integer e,eg,buf(0:49)
#      integer e,eg,buf[0:49]

      nwds = (1 + ldim*(2^ldim))*(wdsizi/4) # group + 2x4 for 2d, 3x8 for 3d

      if ifbswap && ierr == 0 && wdsizi == 8
          byte_reverse8(buf, nwds, ierr)
      elseif ifbswap && ierr == 0 && wdsizi == 4
          byte_reverse(buf, nwds, ierr)
      end
      if (ierr != 0) return end

      if wdsizi == 8
         copyi4(igroup(e), buf[0], 1) #0-1
         if ldim == 3
            COPY(xc(1, e), buf[2], 8) #2 --17
            COPY(yc(1, e), buf[18], 8) #18--33
            COPY(zc(1, e), buf[34], 8) #34--49
         else
            COPY(xc(1, e), buf[2], 4) #2 --9
            COPY(yc(1, e), buf[10], 4) #10--17
          end
      else
         igroup(e) = buf[0]
         if if3d
            copy4r(xc(1, e), buf[1], 8)
            copy4r(yc(1, e), buf[9], 8)
            copy4r(zc(1, e), buf[17], 8)
         else
            copy4r(xc(1, e), buf[1], 4)
            copy4r(yc(1, e), buf[5], 4)
         end
      end

      return nothing
      end
#-----------------------------------------------------------------------
      function buf_to_curve(buf)    # version 1 of binary reader

      include("SIZE")
      include("TOTAL")

#      integer e,eg,f,buf[30]

      if wdsizi == 8
        copyi4(eg, buf[1], 1) #1-2
        e  = gllel(eg)

        copyi4(f, buf[3], 1) #3-4

        COPY( curve(1, f, e), buf[5] , 5) #5--14
        chcopy(ccurve(  f, e), buf[15], 1)#15
      else
        eg = buf[1]
        e  = gllel(eg)
        f  = buf[2]

        copy4r( curve(1, f, e), buf[3], 5)
        chcopy(ccurve(f, e)  , buf[8], 1)
      end

#     write(6,1) eg,e,f,(curve(k,f,e),k=1,5),ccurve(f,e)
#   1 format(2i7,i3,5f10.3,1x,a1,'ccurve')

      return nothing
      end
#-----------------------------------------------------------------------
      function buf_to_bc(cbl, bl, buf)    # version 1 of binary reader

      include("SIZE")
      include("TOTAL")

#      character*3 cbl[6,lelt]
#      REAL         bl[5,6,lelt]

#      integer e,eg,f,buf[30]

      if wdsizi == 8
        copyi4(eg, buf[1], 1) #1-2
        e  = gllel(eg)

        copyi4(f, buf[3], 1) #3-4

        COPY(bl[1,f,e], buf[5], 5) #5--14
        chcopy(cbl[f,e], buf[15], 3)#15-16

        if (nelt >= 1000000 && cbl[f,e] == "P  ")
         copyi4(bl[1,f,e], buf[5], 1) #Integer assign connecting P element
        end

      else
        eg = buf[1]
        e  = gllel(eg)
        f  = buf[2]

        copy4r( bl[1,f,e], buf[3], 5)
        chcopy(cbl[f,e], buf[8], 3)

        if (nelgt >= 1_000_000 && cbl[f,e] == "P  ")
           bl[1,f,e] = buf[3] # Integer assign of connecting periodic element
        end
      end



#     write(6,1) eg,e,f,cbl(f,e),' CBC',nid
#  1  format(2i8,i4,2x,a3,a4,i8)

      return nothing
      end
#-----------------------------------------------------------------------
      function bin_rdmesh(1)(ifbswap)    # version 1 of binary reader

      include("SIZE")
      include("TOTAL")
#      logical ifbswap

#      integer e,eg,buf[55]

      if (nio == 0) println(stdout, "  reading mesh ") end

      nwds = (1 + ldim*(2^ldim))*(wdsizi/4) # group + 2x4 for 2d, 3x8 for 3d
      len  = 4*nwds                          # 4 bytes / wd

      if nwds > 55 || isize > 4
         println(stdout, nid, " Error in bin_rd1_mesh: buf size", nwds, isize)
         exitt
      end

      nekgsync()

      niop = 10
      for k = 1:8
         if (nelgt/niop < 100) @goto L10 end
         niop = niop*10
      end
      @label L10

      ierr  = 0
      ierr2 = 0
      len1  = 4
      for eg = 1:nelgt             # sync NOT needed here

         mid = gllnid(eg)
         e   = gllel(eg)
#ifdef DEBUG
         if (nio == 0 && mod(eg, niop) == 0) println(stdout, eg, " mesh read") end
#endif
         if mid != nid && nid == 0              # read & send

            if ierr == 0
              byte_read(buf, nwds, ierr)
              csend(e, ierr, len1, mid, 0)
              if (ierr == 0) csend(e, buf, len, mid, 0) end
            else
              csend(e, ierr, len1, mid, 0)
            end

         elseif mid == nid && nid != 0          # recv & process

            crecv(e, ierr, len1)
            if ierr == 0
              crecv(e, buf, len)
              buf_to_xyz(buf, e, ifbswap, ierr2)
            end

         elseif mid == nid && nid == 0          # read & process

            if ierr == 0
              byte_read(buf, nwds, ierr)
              buf_to_xyz(buf, e, ifbswap, ierr2)
            end
         end

      end
      ierr = ierr + ierr2
      err_chk(ierr, "Error reading .re2 mesh. Abort. \$")

      return nothing
      end
#-----------------------------------------------------------------------
      function bin_rdcurve(1)(ifbswap) # v. 1 of curve side reader

      include("SIZE")
      include("TOTAL")
#      logical ifbswap

#      integer e,eg,buf[55]
#      REAL rcurve

      nwds = (2 + 1 + 5)*(wdsizi/4) #eg+iside+ccurve+curve(6,:,:) !only 5 in rea
      len  = 4*nwds      # 4 bytes / wd

      if nwds > 55 || isize > 4
         println(stdout, nid, " Error in bin_rd1_curve: buf size", nwds, isize)
         exitt
      end

      nekgsync()

      ierr = 0
      len1 = 4
      if nid == 0  # read & send/process

         if wdsizi == 8
           byte_read(rcurve, 2, ierr)
           if (ifbswap) byte_reverse8(rcurve, 2, ierr) end
           ncurve = rcurve
         else
           byte_read(ncurve, 1, ierr)
           if (ifbswap) byte_reverse(ncurve, 1, ierr) end
         end

         if (ncurve != 0) println(stdout, "  reading curved sides ") end
         for k = 1:ncurve
           if ierr == 0
              byte_read(buf, nwds, ierr)
              if wdsizi == 8
                if (ifbswap) byte_reverse8(buf, nwds-2, ierr) end
                copyi4(eg, buf[1], 1)  #1,2
              else
                if (ifbswap) byte_reverse(buf, nwds-1, ierr) end # last is char
                eg  = buf[1]
              end

              mid = gllnid(eg)
              if mid == 0 && ierr == 0
                 buf_to_curve(buf)
              else
                 if ierr == 0
                   csend(mid, buf, len, mid, 0)
                 else
                    @goto L98
                 end
              end
           else
               @goto L98
           end
         end
         @label L98
         buf_close_out  # notify all procs: no more data

      else               # wait for data from node 0

         ncurve_mx = 12*nelt
         for k = 1:ncurve_mx+1   # +1 to make certain we receive the close-out

            crecv(nid, buf, len)
            if wdsizi == 8
               copyi4(ichk, buf[1], 1)
               if (ichk == 0) @goto L99 end
               buf_to_curve(buf)
            elseif buf[1] == 0
                @goto L99
            else
               buf_to_curve(buf)
            end

         end
         @label L99
         buf_close_out

      end
      err_chk(ierr, "Error reading .re2 curved data. Abort.\$")


      return nothing
      end
#-----------------------------------------------------------------------
      function bin_rdbc(1)(cbl, bl, ifbswap) # v. 1 of bc reader

      include("SIZE")
      include("TOTAL")
#      logical ifbswap

#      character*3 cbl[6,lelt]
#      REAL         bl[5,6,lelt]

#      integer e,eg,buf[55]
#      REAL rbc_max

      nwds = (2 + 1 + 5)*(wdsizi/4)   # eg + iside + cbc + bc(5,:,:)
      len  = 4*nwds      # 4 bytes / wd

      if nwds > 55 || isize > 4
         println(stdout, nid, " Error in bin_rd1_bc: buf size", nwds, isize)
         exitt
      end

      for e = 1:nelt   # fill up cbc w/ default
      for k = 1:6
         cbl[k,e] = "E  "
      end
      end

      nekgsync()
      ierr = 0
      len1 = 4
      if nid == 0  # read & send/process

         if wdsizi == 8
           byte_read(rbc_max, 2, ierr)
           if (ifbswap) byte_reverse8(rbc_max, 2, ierr) end # last is char
           nbc_max = rbc_max
         else
           byte_read(nbc_max, 1, ierr)
           if (ifbswap) byte_reverse(nbc_max, 1, ierr) end # last is char
         end

         if (nbc_max != 0) println(stdout, "  reading bc for ifld", ifield) end
         for k = 1:nbc_max
#           write(6,*) k,' dobc1 ',nbc_max
            if ierr == 0
               byte_read(buf, nwds, ierr)
               if wdsizi == 8
                 if (ifbswap) byte_reverse8(buf, nwds-2, ierr) end
                 copyi4(eg, buf[1], 1) #1&2 of buf
               else
                 if (ifbswap) byte_reverse(buf, nwds-1, ierr) end # last is char
                 eg  = buf[1]
               end
               mid = gllnid(eg)
#              write(6,*) k,' dobc3 ',eg,mid

               if mid == 0 && ierr == 0
                   buf_to_bc(cbl, bl, buf)
               else
#                  write(6,*) mid,' sendbc1 ',eg
                   if ierr == 0
                     csend(mid, buf, len, mid, 0)
                   else
                      @goto L98
                   end
#                  write(6,*) mid,' sendbc2 ',eg
               end
#              write(6,*) k,' dobc2 ',nbc_max,eg
            else
                @goto L98
            end
         end
#        write(6,*) mid,' bclose ',eg,nbc_max
         @label L98
         buf_close_outv # notify all procs: no more data

      else               # wait for data from node 0

         nbc_max = 2*ldim*nelt
         for k = 1:nbc_max+1  # Need one extra !

#           write(6,*) nid,' recvbc1',k
            crecv(nid, buf, len)
#           write(6,*) nid,' recvbc2',k,buf(1)

            if wdsizi == 8
               copyi4(ichk, buf[1], 1)
               if (ichk == 0) @goto L99 end
               buf_to_bc(cbl, bl, buf)
            elseif buf[1] == 0
                 @goto L99
            else
                buf_to_bc(cbl, bl, buf)
            end

         end
         @label L99
         buf_close_outv

      end

      err_chk(ierr, "Error reading boundary data for re2. Abort.\$")

      return nothing
      end
#-----------------------------------------------------------------------
      function buf_close_outv()# this is the stupid O(P) formulation

      include("SIZE")
      include("PARALLEL")
#      integer*4 ZERO
#      REAL      rzero

      len   = wdsizi
      rzero = 0
      ZERO  = 0
#     write(6,*) nid,' bufclose'
      if nid == 0
         for mid = 1:np-1
            if (wdsizi == 8) csend(mid, rzero, len, mid, 0) end
            if (wdsizi == 4) csend(mid, ZERO, len, mid, 0) end
#           write(6,*) mid,' sendclose'
         end
      end

      return nothing
      end
#-----------------------------------------------------------------------
      function buf_close_out()# this is the stupid O(P) formulation

      include("SIZE")
      include("PARALLEL")
#      integer*4 ZERO
#      REAL      rzero

#     len  = 4
      len   = wdsizi
      ZERO = 0
      rzero = 0
      if nid == 0
         for mid = 1:np-1
            if (wdsizi == 8) csend(mid, rzero, len, mid, 0) end
            if (wdsizi == 4) csend(mid, ZERO, len, mid, 0) end
         end
      end

      return nothing
      end
#-----------------------------------------------------------------------
      function read_rehdr(2)(ifbswap) # open file & chk for byteswap

      include("SIZE")
      include("TOTAL")

#      logical ifbswap,if_byte_swap_test

#      integer fnami[33]
#      character*132 fname
      equivalence(fname, fnami)

#      character*132 hdr
#      character*5 version
#      REAL*4      test

#      logical iffound

      ierr = 0

      if nid == 0
         @printf(stdout, "%s%s\n", " Reading ", re2fle)
         izero(fnami, 33)
         m = indx2(re2fle, 132, " ", 1)-1
         chcopy(fname, re2fle, m)

         inquire(file = fname, exist = iffound)
         if (!iffound) ierr = 1 end
      end
      err_chk(ierr, " Cannot find re2 file!\$")

      if nid == 0
         byte_open(fname, ierr)
         if (ierr != 0) @goto L100 end
         byte_read(hdr, 20, ierr)
         if (ierr != 0) @goto L100 end

         READ(hdr, "%5s%9i%3i%9i", version, nelgt, ldimr, nelgv)
#    1    format(a5,i9,i3,i9)

         wdsizi = 4
         if (version == "#v002") wdsizi = 8 end
         if version == "#v003"
           wdsizi = 8
           param(32) = 1
         end

         byte_read(test, 1, ierr)
         if (ierr != 0) @goto L100 end
         ifbswap = if_byte_swap_test(test, ierr)
         if (ierr != 0) @goto L100 end

         byte_close(ierr)
      end

      @label L100
      err_chk(ierr, "Error reading re2 header\$")

      bcast(wdsizi, ISIZE)
      bcast(ifbswap, LSIZE)
      bcast(nelgv  , ISIZE)
      bcast(nelgt  , ISIZE)
      bcast(ldimr  , ISIZE)
      bcast(param(32), WDSIZE)

      if (wdsize == 4 && wdsizi == 8)
         exitti("wdsize=4 & wdsizi(re2)=8 not compatible\$", wdsizi)
      end

      return nothing
      end
#-----------------------------------------------------------------------
